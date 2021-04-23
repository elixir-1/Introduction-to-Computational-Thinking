### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c7c324de-a354-11eb-287e-2d1d6a6ab5e0
begin
    import Pkg
	
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
        Pkg.PackageSpec(name="StatsBase", version="0.33"),
    ])
	
    using Plots, PlutoUI, StatsBase, Statistics
end

# ╔═╡ ef0d7c89-6017-4f42-8a85-db24dafb67e1
PlutoUI.TableOfContents(aside=true)

# ╔═╡ fb7c54c6-79d8-41f9-bb51-1bbb3193588c
md"""
# Modelling a simple Pendulum
"""

# ╔═╡ 4df74a52-c618-4e04-b44b-0d610c3b58a0
md"""
Modelling the motion of a simple pendulum consisting of a stiff, inextensible rod of length l with a bob of mass m. It's position can be determined by the angle $\theta$ that it makes with the central reference line.

Let the gravtitational acceleration be denoted by $g$.

Restorative force on the pendulum is given by,

$$F_R = - m g sin(\theta)$$

This is equal to mass times acceleration a,

$$F = m a$$

$$m a = -m g sin(\theta)$$

$$a = - g sin(\theta)$$

Acceleration is second derivative of displacement s,

$$a = \frac{d^2 s}{dt^2}$$

where

$$s = l \theta$$

Hence

$$a = l \frac{d^2 \theta}{dt^2}$$

Let the velocity be denoted by v then it can be represented as,

$$l \frac{d^2 \theta}{dt^2} = - g sin(\theta)$$

$$\frac{d^2 \theta}{dt^2} = - \frac{g}{l} sin(\theta)$$ 

or 

$$\ddot{\theta} = - \frac{g}{l} sin(\theta)$$
"""

# ╔═╡ d1af15e4-2beb-4a7f-9388-d78523c433a7
md"""
## Euler Method to find pendulum position
"""

# ╔═╡ 3ea4dc17-4868-404f-a575-b97a79badbb1
md"""
Let the initial position of the pendulum be $\theta_0$ and we want to find to the position of the pendulum in terms of the angle formed with the vertical reference line on $k^{th}$ day as $\theta_k$.

The position on the pendulum on day k+1, can be found using the position  

Using the Euler method to discreetize the equation above,

$$\ddot{\theta} = - \frac{g}{l} sin(\theta)$$

Approximate the derivative using a small (but not *too* small) time step $\delta t$:

$$\ddot{\theta} \simeq \frac{\dot{\theta}(t + \delta t) - \dot{\theta}(t)}{\delta t}$$

giving

$$\dot{\theta}(t + \delta t) \simeq \dot{\theta}(t) + \delta t \, \ddot{\theta}$$

and

$$\theta(t + \delta t) \simeq \theta(t) + \delta t \, \dot{\theta}(t)$$
"""

# ╔═╡ 0e051bc8-e67e-4ca5-8bb4-1f238bbe3f7b
md"""
Initial value conditions for the pendulum model are,

Initial position at t=0,

$$\theta(t=0) = \theta_0$$

$$\dot{\theta}(t=0) = 0$$
"""

# ╔═╡ 203db499-963a-4fa3-957f-b905424d96d9
let
	m = 1
	g = 9.8
	l = 10
	
	T = 100
	δt = 1
	
	θ_0 = π/12
	θ̇_0 = 0
	θ̈_0 = - (g / l) * sin(θ_0)
	
	time_range = 0+δt:δt:T
	
	axis_lim = 10
	
	θ = Vector{Float64}()
	append!(θ, θ_0)
	
	x = l * sin(θ_0)
	y = - l * cos(θ_0)
	
	PE = m * g * l * cos(θ_0)
	KE = 0.5 * m * l * l * (θ̇_0 * θ̇_0)
		
	ME = PE + KE
	
	me_t = Vector{Float64}()
	append!(me_t, ME)
	
	# plot([0,x], [0,y],size=(400,300),xlim=(-2*axis_lim,2*axis_lim),ylim=(-1.5*axis_lim,axis_lim),markersize = 8, markershape = :circle,label ="",axis = []);
	
	for t in time_range
		str = string("Time = ", round(t, digits = 2), " sec");
		
		θ_t = θ_0 + (δt * θ̇_0)
		append!(θ, θ_t)
		
		x = l * sin(θ_t)
		y = - l * cos(θ_t)
		
		# plot!([0,x], [0,y],size=(400,300),xlim=(-2*axis_lim,2*axis_lim),ylim=(-1.5*axis_lim,axis_lim),markersize = 8, markershape = :circle,label ="",axis = [], title = str, title_location = :left);
		
		θ̇_0 = θ̇_0 + (δt * θ̈_0)
		θ_0 = θ_t
		θ̈_0 = - (g / l) * sin(θ_0)
		
		PE = m * g * l * cos(θ_0)
		KE = 0.5 * m * l * l * (θ̇_0 * θ̇_0)
		
		ME = PE + KE
		append!(me_t, ME)
		
	end
	
	# plot(θ, alpha=0.5, ms=2, label="")
	plot(me_t, alpha=0.5, ms=2, label="Energy in the simple pendulum")
end

# ╔═╡ 1ac63819-9e52-4e77-b325-109775d4f04d
begin
	T = 17632
	δt = 4
	time_range = 0+δt:δt:T
	x = Vector{Float64}()
	append!(x, 0)
	for t in time_range
		append!(x, t)
	end
	[size(time_range, 1)+1 size(x, 1)]
	0+δt:δt:T
end

# ╔═╡ 3bb8acbf-a700-414f-ac9b-f53847b6fc77
begin
	T1 = 100
	δt1 = 0.1
	time_range1 = 0+δt1:δt1:T1
	x1 = Vector{Float64}()
	for t in time_range1
		append!(x1,t)
		if t>3
			append!(x1, t)
		end
	end
	x1
end

# ╔═╡ Cell order:
# ╠═c7c324de-a354-11eb-287e-2d1d6a6ab5e0
# ╠═ef0d7c89-6017-4f42-8a85-db24dafb67e1
# ╟─fb7c54c6-79d8-41f9-bb51-1bbb3193588c
# ╟─4df74a52-c618-4e04-b44b-0d610c3b58a0
# ╟─d1af15e4-2beb-4a7f-9388-d78523c433a7
# ╠═3ea4dc17-4868-404f-a575-b97a79badbb1
# ╟─0e051bc8-e67e-4ca5-8bb4-1f238bbe3f7b
# ╠═203db499-963a-4fa3-957f-b905424d96d9
# ╟─1ac63819-9e52-4e77-b325-109775d4f04d
# ╟─3bb8acbf-a700-414f-ac9b-f53847b6fc77
