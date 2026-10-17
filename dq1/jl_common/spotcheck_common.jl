# spotcheck_common.jl -- unified ΔQ1 spot check (shared by the four Julia PRs): argument parsing, seeds and uniform sampling, branch rotation,
# tolerance verdict, replay of rejected points, truncation by time budget, JSON output (write .tmp first, then rename). Depends on no non-standard-library package.
module SpotCheck

using Random, TOML, Printf

export Check, chk_close, chk_le, chk_bool, chk_ge, Branch, run_tracks, parse_cli,
       unif, logunif, pick, maxabs, relmaxdiff, pkg_meta, seed_for

#####
##### Verdict primitives
#####

struct Check
    name :: String
    obs  :: Any
    exp  :: Any
    res  :: Float64
    tol  :: Float64
    pass :: Bool
end

_f(x) = x isa Number ? Float64(x) : NaN
_finite(x) = x isa Number ? isfinite(x) : true

"|obs - exp| ≤ atol + rtol·|exp|"
function chk_close(name, obs, exp; rtol = 1e-6, atol = 1e-12)
    r = abs(_f(obs) - _f(exp))
    t = atol + rtol * abs(_f(exp))
    return Check(name, obs, exp, r, t, isfinite(r) && r <= t)
end

"lhs ≤ rhs + tol; res = lhs - rhs"
function chk_le(name, lhs, rhs; tol = 0.0)
    r = _f(lhs) - _f(rhs)
    return Check(name, lhs, rhs, r, Float64(tol), isfinite(r) && r <= tol)
end

"lhs ≥ rhs - tol; res = rhs - lhs"
function chk_ge(name, lhs, rhs; tol = 0.0)
    r = _f(rhs) - _f(lhs)
    return Check(name, lhs, rhs, r, Float64(tol), isfinite(r) && r <= tol)
end

chk_bool(name, pass::Bool; obs = pass, exp = true, res = pass ? 0.0 : 1.0, tol = 0.0) =
    Check(name, obs, exp, Float64(res), Float64(tol), pass)

maxabs(a) = isempty(a) ? 0.0 : Float64(maximum(abs, a))
"max|a-b| / max(max|b|, floor)"
relmaxdiff(a, b; floor = 0.0) = maxabs(Float64.(a) .- Float64.(b)) / max(maxabs(b), floor)

#####
##### Sampling
#####

seed_for(prseq, track) = 20261017 + prseq * 10 + track
unif(rng, a, b) = a + (b - a) * rand(rng)
logunif(rng, a, b) = exp(log(a) + (log(b) - log(a)) * rand(rng))
pick(rng, xs) = xs[rand(rng, 1:length(xs))]

#####
##### Branches
#####

struct Branch
    name    :: String
    sampler :: Function   # rng -> NamedTuple input
    eval    :: Function   # input -> Vector{Check}; an empty vector = nothing to check at this point (recorded in n_empty)
end

#####
##### JSON (ordered objects use Vector{Pair})
#####

_jesc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")
function _jnum(x::Real)
    x isa Bool && return x ? "true" : "false"
    x isa Integer && return string(x)
    y = Float64(x)
    isnan(y) && return "\"NaN\""
    isinf(y) && return y > 0 ? "\"Inf\"" : "\"-Inf\""
    return repr(y)
end
function tojson(io, x)
    if x === nothing
        print(io, "null")
    elseif x isa Real
        print(io, _jnum(x))
    elseif x isa AbstractString || x isa Symbol || x isa Char
        print(io, '"', _jesc(x), '"')
    elseif x isa AbstractVector{<:Pair} || x isa AbstractDict || x isa NamedTuple
        ps = x isa NamedTuple ? collect(pairs(x)) : collect(x)
        print(io, '{')
        for (i, (k, v)) in enumerate(ps)
            i > 1 && print(io, ',')
            print(io, '"', _jesc(k), "\":")
            tojson(io, v)
        end
        print(io, '}')
    elseif x isa AbstractArray || x isa Tuple
        print(io, '[')
        for (i, v) in enumerate(x)
            i > 1 && print(io, ',')
            tojson(io, v)
        end
        print(io, ']')
    elseif x isa DataType || x isa Function
        print(io, '"', _jesc(string(x)), '"')
    else
        print(io, '"', _jesc(string(x)), '"')
    end
end
function write_json_atomic(path, obj)
    mkpath(dirname(abspath(path)))
    tmp = path * ".tmp"
    open(tmp, "w") do io
        tojson(io, obj)
        println(io)
    end
    mv(tmp, path; force = true)
end

#####
##### Version information
#####

function pkg_meta(names)
    proj = Base.active_project()
    man = joinpath(dirname(proj), "Manifest.toml")
    out = Pair{String,Any}["julia" => string(VERSION), "project" => proj]
    isfile(man) || return out
    d = TOML.parsefile(man)
    deps = get(d, "deps", Dict())
    for n in names
        haskey(deps, n) || continue
        e = deps[n][1]
        push!(out, n => Pair{String,Any}[k => string(get(e, k, "")) for k in ("version", "repo-rev", "repo-url", "git-tree-sha1")])
    end
    return out
end

#####
##### Running one track
#####

_inputjson(inp) = inp isa NamedTuple ? Pair{String,Any}[string(k) => (v isa Number || v isa AbstractString || v isa Symbol || v isa Bool ? v : string(v)) for (k, v) in pairs(inp)] : string(inp)
_checkjson(c::Check) = Pair{String,Any}["check" => c.name, "observed" => (c.obs isa Number ? c.obs : string(c.obs)),
                                        "expected" => (c.exp isa Number ? c.exp : string(c.exp)),
                                        "residual" => c.res, "tol" => c.tol]

mutable struct BStat
    n::Int; n_reject::Int; n_fail::Int; n_empty::Int; t::Float64
    checks::Dict{String, Vector{Int}}   # name => [n_evaluated, n_reject, n_nonfinite]
    obsrange::Dict{String, Vector{Float64}}
end
BStat() = BStat(0, 0, 0, 0, 0.0, Dict{String,Vector{Int}}(), Dict{String,Vector{Float64}}())

function _evaluate(b::Branch, inp)
    t0 = time()
    checks = try
        b.eval(inp)
    catch err
        msg = sprint(showerror, err)
        return (:fail, Check[], time() - t0, first(msg, 600))
    end
    dt = time() - t0
    isempty(checks) && return (:empty, checks, dt, "")
    if any(c -> !_finite(c.obs) || !_finite(c.exp) || isnan(c.res), checks)
        return (:fail, checks, dt, "non-finite observable (NaN/Inf): " *
                join([c.name for c in checks if !_finite(c.obs) || !_finite(c.exp) || isnan(c.res)], ","))
    end
    return (all(c -> c.pass, checks) ? :pass : :reject, checks, dt, "")
end

function run_track(; pr, prseq, track, n, out, budget = Inf, branches::Vector{Branch},
                   must = Tuple{String, Any}[], meta = Pair{String,Any}[], max_replay = 20,
                   rotate = length(branches))
    seed = seed_for(prseq, track)
    rng = Xoshiro(seed)
    bidx = Dict(b.name => b for b in branches)
    stats = Dict(b.name => BStat() for b in branches)
    must_stats = Dict(b.name => BStat() for b in branches)
    rejects = Any[]; fails = Any[]
    nrep = 0
    T0 = time()
    function record!(st, b, inp, status, checks, dt, msg, ismust, idx)
        st.n += 1; st.t += dt
        for c in checks
            v = get!(st.checks, c.name, [0, 0, 0])
            v[1] += 1
            (!_finite(c.obs) || !_finite(c.exp) || isnan(c.res)) ? (v[3] += 1) : (!c.pass && (v[2] += 1))
            if c.obs isa Number && isfinite(_f(c.obs))
                r = get!(st.obsrange, c.name, [Inf, -Inf, -Inf])   # obs minimum, obs maximum, res/tol maximum
                r[1] = min(r[1], _f(c.obs)); r[2] = max(r[2], _f(c.obs))
                isfinite(c.res) && (r[3] = max(r[3], c.tol > 0 ? c.res / c.tol : (c.res > 0 ? Inf : 0.0)))
            end
        end
        if status == :empty
            st.n_empty += 1
        elseif status == :fail
            st.n_fail += 1
            length(fails) < 20 && push!(fails, Pair{String,Any}["branch" => b.name, "point" => idx, "must" => ismust,
                                                               "input" => _inputjson(inp), "error" => msg])
        elseif status == :reject
            st.n_reject += 1
            if length(rejects) < 20
                rep = nothing
                if nrep < max_replay
                    nrep += 1
                    s2, c2, _, _ = _evaluate(b, inp)
                    rep = (s2 == :reject) && length(c2) == length(checks) &&
                          all(isequal(_f(a.obs), _f(c.obs)) for (a, c) in zip(c2, checks))
                end
                push!(rejects, Pair{String,Any}["branch" => b.name, "point" => idx, "must" => ismust,
                    "input" => _inputjson(inp), "replay_identical" => rep,
                    "failed_checks" => [_checkjson(c) for c in checks if !c.pass]])
            end
        end
    end
    println("== $(pr) track $(track): seed $(seed), n=$(n), budget=$(budget) s, branches=",
            join([b.name for b in branches], ","), ", must=", length(must)); flush(stdout)
    for (j, (bn, inp)) in enumerate(must)
        b = bidx[bn]
        status, checks, dt, msg = _evaluate(b, inp)
        record!(must_stats[bn], b, inp, status, checks, dt, msg, true, j)
        @printf("  must %d %s %s %.2fs\n", j, bn, status, dt); flush(stdout)
    end
    k = rotate   # rotate only among the first `rotate` branches; later branches have no input to sample and serve only as mandatory points
    ndone = 0; stopped = false; tlast = time()
    for i in 1:n
        if time() - T0 > budget
            stopped = true
            break
        end
        b = branches[mod1(i, k)]
        inp = b.sampler(rng)
        status, checks, dt, msg = _evaluate(b, inp)
        record!(stats[b.name], b, inp, status, checks, dt, msg, false, i)
        ndone = i
        if time() - tlast > 60 || i <= k
            @printf("  [%s t%d] %d/%d  %.0fs  last=%s(%s %.2fs)\n", pr, track, i, n, time() - T0, b.name, status, dt)
            flush(stdout); tlast = time()
        end
    end
    tot = (n = sum(s.n for s in values(stats)), rej = sum(s.n_reject for s in values(stats)),
           fail = sum(s.n_fail for s in values(stats)), empty = sum(s.n_empty for s in values(stats)))
    mtot = (n = sum(s.n for s in values(must_stats)), rej = sum(s.n_reject for s in values(must_stats)),
            fail = sum(s.n_fail for s in values(must_stats)))
    bjson(st) = Pair{String,Any}["n" => st.n, "n_reject" => st.n_reject, "n_fail" => st.n_fail, "n_empty" => st.n_empty,
                                 "seconds" => st.t, "per_point_seconds" => st.n > 0 ? st.t / st.n : nothing,
                                 "checks" => Pair{String,Any}[c => Pair{String,Any}["n" => v[1], "n_reject" => v[2], "n_nonfinite" => v[3],
                                                                   "obs_min" => get(st.obsrange, c, [NaN])[1],
                                                                   "obs_max" => get(st.obsrange, c, [NaN, NaN])[2],
                                                                   "max_res_over_tol" => get(st.obsrange, c, [NaN, NaN, NaN])[3]]
                                                              for (c, v) in sort(collect(st.checks))]]
    obj = Pair{String,Any}[
        "pr" => pr, "track" => track, "seed" => seed,
        "n" => ndone, "n_planned" => n, "rotating_branches" => [b.name for b in branches[1:k]], "stopped_by_budget" => stopped, "budget_seconds" => budget,
        "branches" => Pair{String,Any}[b.name => bjson(stats[b.name]) for b in branches],
        "n_reject_total" => tot.rej, "n_fail_total" => tot.fail, "n_empty_total" => tot.empty,
        "must_points" => Pair{String,Any}["n" => mtot.n, "n_reject" => mtot.rej, "n_fail" => mtot.fail,
            "branches" => Pair{String,Any}[b.name => bjson(must_stats[b.name]) for b in branches if must_stats[b.name].n > 0]],
        "first_rejections" => rejects, "first_failures" => fails,
        "versions" => meta, "elapsed_seconds" => time() - T0,
        "host" => gethostname(), "julia_threads" => Threads.nthreads()]
    write_json_atomic(out, obj)
    @printf("== %s track %d done: n=%d reject=%d fail=%d empty=%d must(n=%d rej=%d fail=%d) %.0fs -> %s\n",
            pr, track, ndone, tot.rej, tot.fail, tot.empty, mtot.n, mtot.rej, mtot.fail, time() - T0, out)
    flush(stdout)
    return obj
end

#####
##### Command line
#####

"Returns (tracks, n, outpaths::Dict{Int,String}, budgets::Dict{Int,Float64})"
function parse_cli(args, prname)
    tracks = Int[]; n = 100; out = nothing; outdir = nothing; budgets = Float64[]; budget = Inf
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--track"; tracks = [parse(Int, args[i+1])]; i += 2
        elseif a == "--tracks"; tracks = parse.(Int, split(args[i+1], ",")); i += 2
        elseif a == "--n"; n = parse(Int, args[i+1]); i += 2
        elseif a == "--out"; out = args[i+1]; i += 2
        elseif a == "--outdir"; outdir = args[i+1]; i += 2
        elseif a == "--budget"; budget = parse(Float64, args[i+1]); i += 2
        elseif a == "--budgets"; budgets = parse.(Float64, split(args[i+1], ",")); i += 2
        else error("unknown argument $a")
        end
    end
    isempty(tracks) && error("--track or --tracks is required")
    outs = Dict{Int,String}()
    for t in tracks
        outs[t] = out !== nothing && length(tracks) == 1 ? out :
                  joinpath(something(outdir, "."), "$(prname)_t$(t).json")
    end
    bud = Dict{Int,Float64}()
    for (j, t) in enumerate(tracks)
        bud[t] = isempty(budgets) ? budget : budgets[j]
    end
    return tracks, n, outs, bud
end

"Run several tracks in one process; if a whole track errors, write a JSON containing only the error and continue with the next track"
function run_tracks(args; pr, prseq, build_track::Function, meta)
    tracks, n, outs, bud = parse_cli(args, pr)
    for t in tracks
        try
            spec = build_track(t)    # (branches=..., must=...)
            run_track(; pr, prseq, track = t, n, out = outs[t], budget = bud[t],
                      branches = spec.branches, must = spec.must, meta,
                      rotate = hasproperty(spec, :rotate) ? spec.rotate : length(spec.branches))
        catch err
            msg = sprint(showerror, err, catch_backtrace())
            println("!! $(pr) track $(t) aborted: ", first(msg, 2000))
            write_json_atomic(outs[t], Pair{String,Any}["pr" => pr, "track" => t, "seed" => seed_for(prseq, t),
                "aborted" => true, "error" => first(msg, 4000), "versions" => meta])
        end
        GC.gc()
    end
end

end # module
