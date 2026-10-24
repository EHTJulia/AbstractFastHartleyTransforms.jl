module AbstractFastHartleyTransformsEnzymeCoreExt

using AbstractFastHartleyTransforms:
    AbstractFastHartleyTransforms,
    DHTPlan,
    is_inplace,
    plan_fht,
    plan_fht!,
    plan_ifht,
    plan_ifht!
using LinearAlgebra: mul!
using EnzymeCore: EnzymeCore, Const, Duplicated
using EnzymeCore.EnzymeRules:
    EnzymeRules, FwdConfigWidth, RevConfigWidth, needs_primal, needs_shadow, AugmentedReturn

# Planning depends on shape and metadata, not on differentiable sample values.
# Create plans outside differentiation when using destructive FFTW planning flags.
for planner in (plan_fht, plan_fht!, plan_ifht, plan_ifht!)
    @eval EnzymeRules.inactive(::typeof($planner), x::Array, args...) = true
end

# Runtime activity may represent a constant array by aliasing its primal and shadow.
# Reject that unsupported case rather than transforming the constant data twice.
function check_active_buffer(config, x)
    if EnzymeRules.runtime_activity(config) && x.val === x.dval
        throw(
            ArgumentError(
                "Hartley rules require an independent derivative buffer; move dynamically constant transforms outside autodiff.",
            ),
        )
    end
    return nothing
end

# Restrict these rules to host arrays. Device compilation uses its own execution path.
function EnzymeRules.forward(
        config::FwdConfigWidth{1},
        ::Const{typeof(*)},
        RT::Type,
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    primal = if needs_primal(config) || is_inplace(p.val)
        p.val * x.val
    else
        nothing
    end
    tangent = if needs_shadow(config) || is_inplace(p.val)
        p.val * x.dval
    else
        nothing
    end
    if needs_primal(config) && needs_shadow(config)
        return Duplicated(primal, tangent)
    elseif needs_primal(config)
        return primal
    elseif needs_shadow(config)
        return tangent
    end
    return nothing
end

function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(*)},
        RT::Type,
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    primal = if needs_primal(config) || is_inplace(p.val)
        p.val * x.val
    else
        nothing
    end
    shadow = is_inplace(p.val) ? x.dval : zero(x.val)
    return AugmentedReturn(
        needs_primal(config) ? primal : nothing,
        needs_shadow(config) ? shadow : nothing,
        shadow,
    )
end

function EnzymeRules.reverse(
        config::RevConfigWidth{1},
        ::Const{typeof(*)},
        RT::Type,
        shadow,
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    # Real Hartley plans and their real-scaled inverses are symmetric.
    # Copy only for mutating plans to preserve the cotangent.
    contribution = p.val * (is_inplace(p.val) ? copy(shadow) : shadow)
    if is_inplace(p.val)
        x.dval .= contribution
    else
        x.dval .+= contribution
        fill!(shadow, zero(eltype(shadow)))
    end
    return (nothing, nothing)
end

function EnzymeRules.forward(
        config::FwdConfigWidth{1},
        ::Const{typeof(mul!)},
        RT::Type,
        y::Duplicated{<:Array},
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    check_active_buffer(config, y)
    tangent = p.val * (is_inplace(p.val) ? copy(x.dval) : x.dval)
    primal = mul!(y.val, p.val, x.val)
    y.dval .= tangent
    if needs_primal(config) && needs_shadow(config)
        return Duplicated(primal, y.dval)
    elseif needs_primal(config)
        return primal
    elseif needs_shadow(config)
        return y.dval
    end
    return nothing
end

function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(mul!)},
        RT::Type,
        y::Duplicated{<:Array},
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    check_active_buffer(config, y)
    primal = mul!(y.val, p.val, x.val)
    return AugmentedReturn(
        needs_primal(config) ? primal : nothing,
        needs_shadow(config) ? y.dval : nothing,
        nothing,
    )
end

function EnzymeRules.reverse(
        config::RevConfigWidth{1},
        ::Const{typeof(mul!)},
        RT::Type,
        tape,
        y::Duplicated{<:Array},
        p::Const{<:DHTPlan},
        x::Duplicated{<:Array},
    )
    contribution = p.val * (is_inplace(p.val) ? copy(y.dval) : y.dval)
    if y.dval === x.dval
        x.dval .= contribution
    else
        # Consume the contribution before clearing a potentially aliased result.
        x.dval .+= contribution
        fill!(y.dval, zero(eltype(y.dval)))
    end
    return (nothing, nothing, nothing)
end

end
