struct BewleyAiyagariStates{F} <: FieldVector{2,F}
    k::F
    e::F
end

@with_kw mutable struct BewleyAiyagariAggregates{F} <: FieldVector{3,F}
    k::F = 1.0
    r::F = 2.0
    w::F = 3.0
end

@with_kw struct BewleyAiyagariParameters{F,AT}
    β::F = 0.99
    θ::F = 1.5
    α::F = 0.35
    δ::F = 0.02
    Z::F = 1.0
    μ::F = 0.0
    ρ::F = 0.99
    σ::F = 0.10
    bounds::NamedTuple = (; k = (0,2500))
    n::NamedTuple = (;k = 10, e = 10)
    indices::NamedTuple = (;k=Index(n.k), e=Index(n.e), e′=Index(n.e))
    ranges::BewleyAiyagariStates{AbstractRange} = BewleyAiyagariStates(range(bounds.k[1],bounds.k[2],n.k),rowenhorst_range(ρ,μ,σ,n.e))
    knots::BewleyAiyagariStates{AT} = BewleyAiyagariStates(collect.(ranges)...)
    interpolation_types::NamedTuple = (; K = (Gridded(Linear()),NoInterp()), Vbar = (Gridded(Linear()),NoInterp()))
    aggregates::BewleyAiyagariAggregates{F} = BewleyAiyagariAggregates()
    transition_matrix::ITensor = itensor(rowenhorst_matrix(ρ,μ,σ,n.e),indices.e,indices.e′)
    update_rule = (new,old) -> 1/6*(new-old)+old
end

function bewley_aiyagari_steady_state(parameters)
    @unpack_BewleyAiyagariParameters parameters

    # preallocate 
    V = zeros(n...)
    Vcopy = similar(V)
    Q = similar(V)
    Qstar = similar(V)
    K = similar(V)
    Kcopy = similar(V)
    Vbar = interpolate!((knots.k,1:n),similar(V),interpolation_types.Vbar)
    
    # model closures
    function v(k_prime,k,e,eidx,Vbar)
        y = aggregates.w*exp(e)+(1-δ)*k+aggregates.r*k
        c =  y-k′
        ((1-β)*c^(1-θ)+β*Vbar(k′,eidx)^(1-θ))^(1/(1-θ))
    end
    vstar = (k,e,eidx,Vbar) -> begin
        @_ golden_search(v(_1,k,e,eidx,_2),bounds.k[1],min(y,bounds.k[2]),Vbar)
    end

    # iteration steps
    valuestep = () -> begin
        coefficients(Vbar) .= itensor(V,indices.k,indices.e)*itensor(Q,indices.e,indices.e′)
        interpolate!(Vbar)
        @_ map!(v,V,(knots.k,knots.e,axes(knots.e,1)),Vbar)
        return norm_copy!(Vcopy,V)
    end
    
    policystep = () -> begin
        @_ map!(vstar,(V,K),(knots.k,knots.e,axes(knots.e,1)),Vbar)
        return nothing 
    end
    
    ergodic_distribution_step = () -> begin
        copyto!(Kcopy,K)
        K_itp = interpolate!((nodes.k,1:n),Kcopy,interpolation_types.K)
        op = FunctionOperator(similar(Q)) do (v,u,p,t)
            transition!(Qstar,Q,K_itp)
            Q .= itensor(Qstar,indices.k,indices.e) * transition_matrix
        end
        Q .= solve(LinearProblem(op,Q))
        return nothing
    end
    
    aggregate_step = () -> begin
        @unpack_BewleyAiyagariAggregates aggregates
        knew = itensor(Q,indices.k,indices.e)*itensor(K,indices.k,indices.e)
        k = update_rule(knew,k)
        r = α*z*k
        w = (1-α)*z*k
        @pack_BewleyAiyagariAggregates! aggregates
        return nothing
    end
    
end

function bewley_aiyagari_shooting_method(parameters,type)
    @unpack_BewleyAiyagariParameters parameters

    # preallocate
    
    T = itensor(rowenhorst_matrix(ρ,μ,σ,n.e),indices.e,indices.e_prime)
    V = zeros(n...)
    Q = similar(V)
    Qstar = similar(V)
    K = similar(V)
    Vcopy = similar(V)
    Vbar = interpolate!(Tuple(ranges...),similar(V),interpolation_types.Vbar)

    # model closures
    function v(k_prime,k,e,Vbar)
        y = aggregates.w*exp(e)+(1-δ)*k+aggregates.r*k
        c =  y-k′
        ((1-β)*c^(1-θ)+β*Vbar(k′,e)^(1-θ))^(1/(1-θ))
    end
    vstar = (k,e,Vbar) -> begin
        @_ golden_search(v(_1,k,e,_2),bounds.kmin,min(y,bounds.kmax),Vbar)
    end

    # iteration steps
    valuestep = () -> begin
        Vbar .= itensor(V,k,e)*Q
        @_ map!(v,V,(K,knots...),Vbar)
        return norm_copy!(Vcopy,V)
    end
    
    policystep = () -> begin
        @_ map!(vstar,(K,V),Tuple(knots),Vbar)
        return nothing 
    end
    
    ergodic_distribution_step = () -> begin
        K = interpolate!(Tuple(ranges...),similar(V),interpolation_types.K)
        op = FunctionOperator(similar(Q)) do (v,u,p,t)
            transition!(Qstar,Q,K)
            Q .= itensor(Qstar,indices.k,indices.e) * T
        end
        Q .= solve(LinearProblem(op,Q))
    end
    
    aggregatesstep = () -> begin
        @unpack_BewleyAiyagariAggregates aggregates
        knew = itensor(Q,indices.k,indices.e)*itensor(K,indices.k,indices.e)
        k += 1/6*(knew-k)
        r = α*z*k
        w = (1-α)*z*k
        @pack_BewleyAiyagariAggregates! aggregates
    end
    
    shooting_backward_step = () -> begin
        
    end
    
    shooting_forward_step = () -> begin
        
    end
    
    krussel_smith = () -> begin
        
    end
    
    naive_euler_equation_step = () -> begin
        
    end

    euler_equation_step = () -> begin
        
    end
    
end

# bisection and broyden
function bewley_aiyagari_krussel_smith(parameters,type;T=500)
    @unpack_BewleyAiyagariParameters parameters

    V = zeros(n...,T)
    Q = similar(V)
    Qstar = similar(V)
    K = similar(V)
    Vcopy = similar(V)
    Vbar = interpolate!(Tuple(ranges...),similar(V),interpolation_types.Vbar)

    record_backward = () -> begin
    end
    
    record_forward = () -> begin
    end
    
    simlulate_step = () -> begin
        
    end
    
end
