### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 55fcfaca-03b9-4d82-822c-5ddf22bc350e
using CSV, DataFrames, Query

# ╔═╡ 61c1412d-7773-4a26-ba62-c3e91ff5d908
using Statistics, StatsBase, HypothesisTests

# ╔═╡ 859b2e3f-7c97-4d3f-b3a0-cb19f5da6a8a
using Distributions

# ╔═╡ c7fc5e4b-c1bc-4e65-a126-34b9dacce8c3
using Downloads

# ╔═╡ 33a32461-9f78-48e0-8e33-9cf767c18d9b
using PlutoUI, PlutoTest

# ╔═╡ 980de46b-4839-4fd7-b4fe-4ebe1b19bf43
begin
	using Plots
	using LaTeXStrings
	using Printf
	#using PlotlyJS
	plotly()
end

# ╔═╡ 88d46a9d-adda-4750-87e8-58c1b50176cd
md"# Transit Duration Analysis of Kepler Sample"

# ╔═╡ cb9f36fe-8bae-4e54-ae9c-d2fdc466c425
md"## Full Sample"

# ╔═╡ dc9d19bc-cdcb-4a19-bbda-71203b0832a3
md"> ⚠ Some thing is fishy here.  Need to update stellar densities in input file!"

# ╔═╡ f8fa6349-04be-4f6c-9465-85208dfc198d
md"## Checking thresholds choices"

# ╔═╡ 05572fad-3d61-4a99-8511-6915b0cda15a
md"### Host Star Temperature"

# ╔═╡ 08aeb059-fe92-455d-a726-577329d68be7
md"""
Minimum Teff: $(@bind min_teff NumberField(3200.0:50:9200, default=4000.0))
Maximum Teff: $(@bind max_teff NumberField(3900.0:50:9400, default=6600.0))
Number Teff Bins: $(@bind num_teff_bins NumberField(2:10, default=6))
"""

# ╔═╡ 77de3d40-428f-463e-88a7-3eb209b52747
md"### Minimum SNR"

# ╔═╡ 1e682cd4-2ac0-43d2-98df-a054520fd5c5
md"""
Alternative Minium SNR: $(@bind alt_min_SNR NumberField(7.0:0.5:30, default= 10.0))
"""

# ╔═╡ 33d289e5-a205-4a82-a0e6-0b09d50a1f1f
md"### Maximum Impact Parameter"

# ╔═╡ 4276feec-2e6f-4a89-89be-7a152a3833dc
md"""
Alternative Maximum b: $(@bind alt_max_b NumberField(0.0 : 0.01 :1, default=0.9))
"""

# ╔═╡ dac4d650-c8b7-4a24-8420-a0ebc56dccca
md"## Fitting Population to analytic models"

# ╔═╡ da90b894-32da-4187-b1a1-6548bf5b5b12
md"""
> ⚠ Need to update stellar densities in input file!
>
> Wait until fix issue with stellar densities before fitting a model!
"""

# ╔═╡ 561c5fee-e103-4de8-ab47-e11f32a2a9d0
md"""
Main e scale: $(@bind e_for_rayleigh NumberField(0:0.005:1.0, default=0.005))
Include noise? $(@bind sim_noise CheckBox(default=true))
Redraw: $(@bind redraw Button("Redraw"))

Lo e scale: $(@bind e_lo_for_rayleigh NumberField(0:0.001:1.0, default=0.002))
Hi e scale: $(@bind e_hi_for_rayleigh NumberField(0:0.02:1.0, default=0.2))
Fraction in hi e componet: $(@bind frac_hi_e NumberField(0:0.02:1.0, default=0.1))
"""

# ╔═╡ b5d1da7a-c483-42f7-9b64-e31bedf768ba
md"""
Minium SNR: $(@bind min_SNR NumberField(7.0:0.5:30, default= 15.0))
"""

# ╔═╡ d3b2e526-8838-460c-ad02-a6bde5c5ef4f
md"Max b: $(@bind max_b NumberField(0:0.01:1.0, default=0.0) )   (0.0=> max b = 1-Rp/R⋆)"

# ╔═╡ 7fe5f5af-2085-49ab-9a6c-990548ccb5f3


# ╔═╡ 5af5ef4f-a8a4-48b8-8111-c30765cf0ac0
md"## Comparing Subsets"

# ╔═╡ ab590e37-4c5b-4609-876a-c09ee92985c8
md"### Multiplicity"

# ╔═╡ 9f930a87-597d-4997-9819-6a41acdd53bb
md"### Planet Size"

# ╔═╡ 616bb0a9-f6da-458a-bc25-ad3dbcd0d77b
md"### Orbital Period"

# ╔═╡ 4cf1d6c9-ff78-4056-be09-2279ccb3b54c
num_P_bins = 4

# ╔═╡ c859043d-d7ba-47cd-9035-319cdc200720
md"### Host Star Temperature"

# ╔═╡ 4ee2936e-d41f-41f1-b307-60d606629b5b
md"## Null Hypothesis Tests"

# ╔═╡ 98a95e8f-bc1d-4a2a-9665-0e661673b025
md"### Stellar Temperature"

# ╔═╡ ff82c99c-b489-43a8-abc4-ceaf694aabf8
md"### Singles vs Multis"

# ╔═╡ 639dfbeb-60e5-4f59-a3dc-81a0884f0ac4
md"### Planet Size"

# ╔═╡ ae6f7834-f939-4a09-9bd6-bf105c6b05c5
md"### Orbital Period"

# ╔═╡ d498d791-9763-47e5-b5cc-3a41d930dbe7
md"## Geneate Simulated τ Distributions"

# ╔═╡ 893a58a5-f88b-4a59-bf21-c66b7e9404d4
repeat([1,2,3],3)[rand(Bool,9)]

# ╔═╡ 46749afb-7cdb-4956-b835-8623b7cd856a
function calc_taus_assuming_e_rayleigh(df::DataFrame, σ::Real; obs_noise::Bool = false)
	max_b::Union{Float64,Int64}
	n = size(df,1)
	b = rand(n).* (max_b>0 ? max_b : 1.0 .- df.RpRs )
	e = rand.(truncated.(Rayleigh(σ),0.0,1.0.-3.0.*0.005.*df.rstar))
	ω = 2π.*rand(n)
	tau = sqrt.((1.0.+b).*(1.0.-b)).*sqrt.((1.0 .+ e).*(1.0 .- e))./(1.0 .+ e.*cos.((π/2).-ω))
	if obs_noise 
		sigma_rho_star = copy(df.rhostar_ep)
		mask = rand(Bool,n)
		sigma_rho_star[mask] .= df.rhostar_em[mask]
		rho_star = rand.(TruncatedNormal.(df.rhostar,abs.(sigma_rho_star),0.0, Inf))
		P = df.TTVperiod .+ df.TTVperiod_e .* randn(n)
		tau .*= (P./df.TTVperiod.*df.rhostar./rho_star).^(1//3)
	end
	weights = ((1.0 .+ e).*(1.0 .- e))./(1.0 .+ e.*cos.((π/2).-ω))
	ecdf(tau,weights=weights)
end

# ╔═╡ bc6d7cfc-b8da-43ac-84a5-73e323323701
#=
# Commented out to eliminate possible source of errors while debugging stellar density issue.  Can try again once that's fixed. 
function calc_taus_assuming_e_rayleigh(df::DataFrame, σ::Real; obs_noise::Bool = false, num_catalogs::Integer = 1)
	max_b::Union{Float64,Int64}
	n = size(df,1)*num_catalogs
	b = rand(n).* (max_b>0 ? max_b : 1.0 .- repeat(df.RpRs,num_catalogs) )
	e = rand.(truncated.(Rayleigh(σ),0.0,1.0.-3.0.*0.005.*repeat(df.rstar,num_catalogs)))
	ω = 2π.*rand(n)
	tau = sqrt.((1.0.+b).*(1.0.-b)).*sqrt.((1.0 .+ e).*(1.0 .- e))./(1.0 .+ e.*cos.((π/2).-ω))
	if obs_noise 
		sigma_rho_star = repeat(df.rhostar_ep,num_catalogs)
		mask = rand(Bool,n*num_catalogs)
		sigma_rho_star[mask] .= repeat(df.rhostar_em,num_catalogs)[mask]
		rho_star = rand.(TruncatedNormal.(repeat(df.rhostar,num_catalogs),abs.(sigma_rho_star),0.0, Inf))
		P = repeat(df.TTVperiod,num_catalogs) .+ repeat(df.TTVperiod_e,num_catalogs) .* randn(n*num_catalogs)
		tau .*= (P./repeat(df.TTVperiod,num_catalogs).*repeat(df.rhostar,num_catalogs)./rho_star).^(1//3)
	end
	weights = ((1.0 .+ e).*(1.0 .- e))./(1.0 .+ e.*cos.((π/2).-ω))
	ecdf(tau,weights=weights)
end
=#

# ╔═╡ cd152b2a-1690-42d5-b76a-cd71b636aa01
function calc_taus_assuming_e_rayleigh(df::DataFrame, σ::AbstractVector{T}, weights::AbstractVector{T}; obs_noise::Bool = false) where T<:Real
	@assert size(σ,1) ==  size(weights,1)
	@assert sum(weights) ≈ 1.0
	n_pops = size(σ,1)
	ecdfs = Vector{ECDF{typeof(σ)}}(undef,n_pops)
	for i in 1:n_pops
		ecdfs[i] = calc_taus_assuming_e_rayleigh(df,σ[i],obs_noise=obs_noise)
	end
	v = mapreduce(x->x.sorted_values,vcat,ecdfs)
	w = mapreduce(i->ecdfs[i].weights.values .* weights[i],vcat,1:size(ecdfs,1))
	ecdf_combo = ecdf(v,weights=w)
	return ecdf_combo
end


# ╔═╡ 7eb040e9-300b-41d2-9171-9f48c9b98758
md"# Setup"

# ╔═╡ 49c7d182-fd27-4d51-a7d2-c2153beda6af
TableOfContents(aside=true)

# ╔═╡ 2704c820-d307-4acc-869a-e8670232ec38
md"### Read data"

# ╔═╡ 35d1e842-496b-11ec-2e35-17c2d536f22f
begin
	path = joinpath(homedir(),"Documents")
	if !isdir(path) path = homedir() end
	datafile = "koiprops_20211110.csv"
	if ! (filesize(joinpath(path,datafile)) >0)
		download("https://www.dropbox.com/s/yberp4e3pbggvk1/koiprops_20211110.csv?dl=1", output=joinpath(path,datafile))
	end
	@test filesize(joinpath(path,datafile)) >0
end

# ╔═╡ cc668221-1b44-4fcb-83bf-d7f851c8bcbd
df_all = CSV.read(joinpath(path,datafile), DataFrame);

# ╔═╡ 585c1c7a-56bd-49b7-b075-b2232072c893
md"### Filter planet table"

# ╔═╡ 591cebe7-9aa4-4617-9aee-c10c038a7002
md"## Diagnosing problem with rhostar"

# ╔═╡ b3447847-e7f2-4844-be44-01669a0f0e74
names(df_all)

# ╔═╡ 8f46b232-4978-4ea5-bade-d8446e37b8ff
names(df_all)[14]

# ╔═╡ 2f44dc3d-00f9-4fb9-b4af-e15eead627c2
names(df_all)[65]

# ╔═╡ 52ca81e0-fa57-4957-a91d-4b6919c3366e
rhosol = 1.989e33/(4π*6.9634e10^3/3)

# ╔═╡ 5d9c5e6e-d42e-4bf0-99a0-fc89612c15d8
function calc_tcirc(d)
	@assert haskey(d,:TTVperiod)
	@assert haskey(d,:rhostar)
	rhosol = 1.989e33/(4π/3*6.9634e10^3)
	days_in_year = 365.2425
	tdur_circ_diameter_earth = 12.9722 
	tdur_circ = tdur_circ_diameter_earth .* (d.TTVperiod./days_in_year).^(1//3).* (rhosol./d.rhostar).^(1//3)
end

# ╔═╡ 3a0a2bd4-32f4-47d7-a34d-b22880c91395
begin
	df_pl_all = df_all |> 
		@filter(_.StatusFlag[1] == 'P') |>
   		@filter( _.stellar_source ∈ (1,2) ) |>
		@filter( _.Per_lc > 0.0  ) |>
		@filter( _.TTVperiod > 0.0  ) |>
   		DataFrame
	df_pl_all.koi_base = floor.(Int64,df_pl_all.KOI)
	df_systems = df_pl_all |> @groupby(_.koi_base) |> 
		@map({koi_base=key(_),num_pl=length(_), teff=_.teff }) |> DataFrame
	df_pl_all.tdur_circ = calc_tcirc.(eachrow(df_pl_all))
	df_pl_all = innerjoin(df_pl_all,df_systems,on=:koi_base,makeunique=true)
end;

# ╔═╡ 98e17af2-f0b4-4a2a-a9ba-2040da85af88
begin
	if 7<=min_SNR<=30
		df_pl_no_teffcut = df_pl_all |> 
			@filter( _.SNR >= min_SNR  ) |>
			@filter( _.b < (max_b>0 ? max_b : 1.0 - _.RpRs )) |> 
			DataFrame
		make_plots = true
	else
		make_plots = false
		md"Update `min_SNR` before recalculating..."
	end
end;

# ╔═╡ e4581234-6052-4ec4-b6d1-94d3c3c10a28
begin
	if (3200 <= min_teff <= 9300) && (3200 <= max_teff <= 9300)
		df_pl = df_pl_no_teffcut |>	@filter( min_teff <= _.teff <= max_teff ) |> DataFrame
		make_kraft_plot = true
	else
		make_kraft_plot = false
		md"Update `max_teff` before recalculating"
	end
end;

# ╔═╡ 69f7e8be-3541-4fa7-b74b-7ebfc4a0fefc
begin
	h_pl = fit(Histogram, df_pl.tduravg, nbins=200)
	#histogram(df_pl.tduravg, nbins=200)
	plot(h_pl,color=1, label=:none)
	plot!(h_pl.edges,maximum(h_pl.weights).*ecdf(df_pl.tduravg).(h_pl.edges), label=:none, color=1)
	xlabel!("Transit Duration (hr)")
	ylabel!("Count")
end

# ╔═╡ 92ef31db-67ae-4d4c-8c91-69734ee64d7a
begin
	tau_all = df_pl.tduravg./df_pl.tdur_circ
	ntdur_all_cum = ecdf(tau_all)
	x_plt = range(minimum(ntdur_all_cum), stop=maximum(ntdur_all_cum), length=200)
end;

# ╔═╡ eadfde5c-9124-4c9f-87b6-a5c79d50d772
begin
	h_tau_norm = fit(Histogram,tau_all, nbins=200)
	plot(h_tau_norm, label=:none, color=1)
	plot!(h_tau_norm.edges, maximum(h_tau_norm.weights).*ntdur_all_cum.(h_tau_norm.edges), label=:none, color=1)
	xlabel!("Normalized Transit Duration")
	ylabel!("Count")
end

# ╔═╡ 4c77a70d-c003-48eb-a683-359446950878
begin
	df_pl_alt_snr = df_pl_all |> 
			@filter( _.SNR >= alt_min_SNR  ) |>
			@filter( _.b < (max_b>0 ? max_b : 1.0 - _.RpRs )) |> 
			DataFrame
	plt_alt_snr = plot(legend=:bottomright)
	plot!(plt_alt_snr,h_tau_norm.edges,ntdur_all_cum.(h_tau_norm.edges), label="min SNR="*string(min_SNR) )
	plot!(plt_alt_snr,h_tau_norm.edges,ecdf(df_pl_alt_snr.tduravg./df_pl_alt_snr.tdur_circ).(h_tau_norm.edges), label="min SNR="*string(alt_min_SNR) )
	xlabel!(plt_alt_snr,"Normalized Transit Duration")
	ylabel!(plt_alt_snr,"Cumulative")
	xlims!(plt_alt_snr,0,2)
end

# ╔═╡ 8dc05082-c2bd-4120-b7f0-499287953056
begin
	df_pl_alt_maxb = df_pl_all |> 
			@filter( _.SNR >= min_SNR  ) |>
			@filter( _.b < (alt_max_b>0 ? alt_max_b : 1.0 - _.RpRs )) |> 
			DataFrame
	plt_alt_maxb = plot(legend=:bottomright)
	plot!(plt_alt_maxb,h_tau_norm.edges,ntdur_all_cum.(h_tau_norm.edges), label="max b="*string(max_b) )
	plot!(plt_alt_maxb,h_tau_norm.edges,ecdf(df_pl_alt_maxb.tduravg./df_pl_alt_maxb.tdur_circ).(h_tau_norm.edges), label="max b SNR="*string(alt_max_b) )
	xlabel!(plt_alt_maxb,"Normalized Transit Duration")
	ylabel!(plt_alt_maxb,"Cumulative")
	xlims!(plt_alt_maxb,0,2)
end

# ╔═╡ 56effba4-4537-4f1b-8225-9683668ca37a
if make_plots
	redraw
	plt_rayleighs = plot(xlabel="τ", ylabel="Cumulative", legend=:topleft)
	plot!(plt_rayleighs,(x_plt),ntdur_all_cum.(x_plt), label="Kepler", color=:lightgrey) 
	plot!(plt_rayleighs,(x_plt.*rhosol^(1//3)),ntdur_all_cum.(x_plt), label="Scaled Kepler", color=:black) 
	plot!(plt_rayleighs,(x_plt),calc_taus_assuming_e_rayleigh(df_pl,e_for_rayleigh,obs_noise=sim_noise).(x_plt), label="R($e_for_rayleigh)", color=:cyan)
	#plot!(plt_rayleighs,(x_plt),calc_taus_assuming_e_rayleigh(df_pl,e_lo_for_rayleigh,obs_noise=sim_noise).(x_plt), label="R($e_lo_for_rayleigh)", color=:red)
	#plot!(plt_rayleighs,(x_plt),calc_taus_assuming_e_rayleigh(df_pl,e_hi_for_rayleigh,obs_noise=sim_noise).(x_plt), label="R($e_hi_for_rayleigh)", color=:blue)
	plot!(plt_rayleighs,(x_plt),calc_taus_assuming_e_rayleigh(df_pl,[e_lo_for_rayleigh,e_hi_for_rayleigh],[1-frac_hi_e,frac_hi_e],obs_noise=sim_noise).(x_plt), label="Mixture Model", color=:red)
	xlims!(plt_rayleighs,0.25,1.5)
end

# ╔═╡ 5eccb812-d7bb-4f9f-a7a0-6892d8eb46ac
begin
	P_limits = quantile(df_pl.TTVperiod,range(0.0,stop=1.0,length=num_P_bins+1))
	df_P_bins = Vector{DataFrame}(undef,num_P_bins)
	for i in 1:num_P_bins
		df_P_bins[i] = df_pl |> @filter( P_limits[i] <= _.TTVperiod <= P_limits[i+1] ) |> DataFrame
	end
	df_P_bins
end;

# ╔═╡ 3dfa1096-7625-409c-a3a3-24f3907d4ac5
begin
	plt_P_bins = plot(legend=:bottomright,palette=:Paired_10)
	for i in 1:length(df_P_bins)
		plot!(plt_P_bins,h_tau_norm.edges,ecdf(df_P_bins[i].tduravg./df_P_bins[i].tdur_circ).(h_tau_norm.edges), label=string(round(Int64,P_limits[i]*10)/10) *" - " * string(round(Int64,P_limits[i+1]*10)/10) * " d", color=i )
	end
	xlabel!(plt_P_bins,"Normalized Transit Duration")
	ylabel!(plt_P_bins,"Cumulative")
	xlims!(plt_P_bins,0,2)
end

# ╔═╡ 3013ec6b-3139-4fdd-aa69-fe4b45dfd9f8
begin
	ad_P_bins = KSampleADTest(Tuple(map(tb->df_P_bins[tb].tduravg./df_P_bins[tb].tdur_circ,1:length(df_P_bins)))... )
end

# ╔═╡ f8cc2102-eeb2-4b4c-a8d0-ec5c206a7f42
md"k-sample Anderson-Darling test:  p-value =  $(pvalue(ad_P_bins))"

# ╔═╡ 1a7dc30e-a793-435c-9c79-c96ad455bb11
P_limits

# ╔═╡ 74fd8933-74f4-4364-a2f8-1465740987b9
begin
	num_Rp_bins = 8
	#rp_limits = quantile(df_pl.radius,range(0.0,stop=1.0,length=num_Rp_bins+1))
	rp_limits = [0.0, 1.0, 1.6, 2.0, 2.5, 3.0, 6.0, 20.0]; num_Rp_bins = length(rp_limits)-1
	df_Rp_bins = Vector{DataFrame}(undef,num_Rp_bins)
	for i in 1:num_Rp_bins
		df_Rp_bins[i] = df_pl |> @filter( rp_limits[i] <= _.radius <= rp_limits[i+1] ) |> DataFrame
	end
	df_Rp_bins
end;

# ╔═╡ 95dfd011-b505-4c94-93f8-9968e5eb6cab
begin
	plt_Rp_bins = plot(legend=:bottomright,palette=:Paired_7)
	for i in 1:length(df_Rp_bins)
		plot!(plt_Rp_bins,h_tau_norm.edges,ecdf(df_Rp_bins[i].tduravg./df_Rp_bins[i].tdur_circ).(h_tau_norm.edges), label=string(round(Int64,rp_limits[i]*10)/10) *" - " * string(round(Int64,rp_limits[i+1]*10)/10) * " R_⊕", color=i )
	end
	xlabel!(plt_Rp_bins,"Normalized Transit Duration")
	ylabel!(plt_Rp_bins,"Cumulative")
	xlims!(plt_Rp_bins,0,2)
end

# ╔═╡ bda80d86-3430-45ef-9006-d6b77b1758d8
begin
	ad_Rp_bins = KSampleADTest(Tuple(map(tb->df_Rp_bins[tb].tduravg./df_Rp_bins[tb].tdur_circ,1:length(df_Rp_bins)))... )
end

# ╔═╡ df6285f5-9e8a-4dee-bbd3-64f1915ed5ba
md"k-sample Anderson-Darling test:  p-value =  $(pvalue(ad_Rp_bins))"

# ╔═╡ 0ac069bb-3e86-4648-964e-69c640538187
begin
scatter(df_pl.tdur_circ,df_pl.tduravg,ms=1.0)
plot!(df_pl.tdur_circ,df_pl.tdur_circ)
end

# ╔═╡ ed49e14f-4b9c-42d2-8e55-2b3091fe4c14
let
	plt =plot(xscale=:log10)
	scatter!(plt,df_pl.TTVperiod,log10.(df_pl.tduravg./df_pl.TTVperiod.^(1//3)),ms=1.0)
	scatter!(plt,df_pl.TTVperiod,log.(df_pl.tdur_circ./df_pl.TTVperiod.^(1//3)), ms=1.4)
end

# ╔═╡ 9179cdfd-1f7a-476f-a535-bae6a2fc192f
scatter(df_pl.rhostar.*rhosol./df_pl.rhostar_model,yerr=df_pl.rhostar_ep.*rhosol./df_pl.rhostar_model,ylims=(0,2),ms=1.0)

# ╔═╡ 0c3d9129-c9be-4dc4-a1a0-1b8ec39ba16c
let
	plt =plot(xscale=:log10)
	histogram!(plt,df_pl.tdur,alpha=0.5, nbins=400)
	histogram!(plt,df_pl.tdur_circ./rhosol, alpha=0.5, nbins=400)
	xlims!(0.2,10)
end

# ╔═╡ cc173a8d-09a4-4ff0-a695-6eba16ade3b1
begin
	teff_round_unit = 100
	teff_limits = round.(quantile(df_pl.teff,range(0.0,stop=1.0,length=num_teff_bins+1))/teff_round_unit,sigdigits=2)*teff_round_unit
	# If prefer to hard code a list of Teff boundaries 
	#teff_limits = [3900.0, 4400, 5300, 6000, 6200,6500, 6600, 6800]; num_teff_bins = length(teff_limits)-1
#	teff_limits = [3264.0, 3900, 3989.0, 4364, 4739, 5114, 5488, 5863, 6238, 6612,7000]; num_teff_bins = length(teff_limits)-1
	df_teff_bins = Vector{DataFrame}(undef,num_teff_bins)
	for i in 1:num_teff_bins
		df_teff_bins[i] = df_pl_no_teffcut |> @filter( teff_limits[i] <= _.teff <= teff_limits[i+1] ) |> DataFrame
	end
	df_teff_bins
end;

# ╔═╡ 3d97e21f-9785-46d8-aa61-87a69fda8b19
begin
	plt_teff_bins = plot(legend=:bottomright,palette=:Paired_10)
	for i in 1:length(df_teff_bins)
		plot!(plt_teff_bins,h_tau_norm.edges,ecdf(df_teff_bins[i].tduravg./df_teff_bins[i].tdur_circ).(h_tau_norm.edges), label=string(round(Int64,teff_limits[i])) *" - " * string(round(Int64,teff_limits[i+1])), color=i )
	end
	xlabel!(plt_teff_bins,"Normalized Transit Duration")
	ylabel!(plt_teff_bins,"Cumulative")
	xlims!(plt_teff_bins,0,2)
end

# ╔═╡ 5cf370f9-875f-4298-bcb2-7481473b2a41
begin
	ad_teff_bins = KSampleADTest(Tuple(map(tb->df_teff_bins[tb].tduravg./df_teff_bins[tb].tdur_circ,1:length(df_teff_bins)))... )
end

# ╔═╡ d80424df-bcfa-4bf1-805e-954fd4c9c25c
md"k-sample Anderson-Darling test:  p-value =  $(round(pvalue(ad_teff_bins),sigdigits=2))"

# ╔═╡ 860b0cdc-67d3-47ce-856f-10f172a97345
begin
	tb = 1
	ks_teff_1_2 = KSampleADTest(df_teff_bins[tb].tduravg./df_teff_bins[tb].tdur_circ, df_teff_bins[tb+1].tduravg./df_teff_bins[tb+1].tdur_circ )
end

# ╔═╡ ea72041d-6bdf-4263-9a78-884edaf2c0af
teff_limits

# ╔═╡ ac1dcf04-19e5-4221-891f-99747b9cf5f4
teff_limits

# ╔═╡ 9dbf4174-17c6-421f-ba9e-953221fcd58f
md"### Create Subsets of planet table"

# ╔═╡ b97166df-bf9a-4208-ac69-6475e6e1cdfa
begin
	df_pl_single = df_pl |> @filter(_.num_pl==1) |> DataFrame;
	df_pl_double = df_pl |> @filter(_.num_pl==2) |> DataFrame;
	df_pl_multi = df_pl |> @filter(_.num_pl>=2) |> DataFrame;
	df_pl_himulti = df_pl |> @filter(_.num_pl>=3) |> DataFrame;
end;

# ╔═╡ a6074ab0-eacd-427b-971a-20ee9af02abb
if make_plots
	plt_multis = plot(xlabel="Normalized Transit Duration", ylabel="Cumulative", legend=:topleft)
	plot!(plt_multis,x_plt,ntdur_all_cum.(x_plt), label="All") 
		
	ntdur_single_cum = ecdf(df_pl_single.tduravg./df_pl_single.tdur_circ)
	plot!(plt_multis,x_plt, ntdur_single_cum.(x_plt), label="Singles")
	ntdur_double_cum = ecdf(df_pl_double.tduravg./df_pl_double.tdur_circ)
	plot!(plt_multis,x_plt, ntdur_double_cum.(x_plt), label="Doubles")
	ntdur_himulti_cum = ecdf(df_pl_himulti.tduravg./df_pl_himulti.tdur_circ)
	plot!(plt_multis,x_plt, ntdur_himulti_cum.(x_plt), label="3+")
	#title!("Grouped by Multiplicity")
	xlims!(0,2)
end

# ╔═╡ 18dbba04-1be6-4a59-9d11-c2c95982897b
ks_single_multi = ApproximateTwoSampleKSTest(df_pl_single.tduravg./df_pl_single.tdur_circ, df_pl_multi.tduravg./df_pl_multi.tdur_circ )

# ╔═╡ 15f7a4d8-d827-441a-b599-902b57f19821
ad_single_multi = KSampleADTest(df_pl_single.tduravg./df_pl_single.tdur_circ, df_pl_multi.tduravg./df_pl_multi.tdur_circ )

# ╔═╡ 9a076cb1-f5fa-45fd-a624-071c2677691d
ks_single_double = ApproximateTwoSampleKSTest(df_pl_single.tduravg./df_pl_single.tdur_circ, df_pl_double.tduravg./df_pl_double.tdur_circ )

# ╔═╡ 2cd616cd-778e-490a-9826-343cf3b42c9d
ad_single_double = KSampleADTest(df_pl_single.tduravg./df_pl_single.tdur_circ, df_pl_double.tduravg./df_pl_double.tdur_circ )

# ╔═╡ ad551ae3-e000-4259-98ea-4cad99f6b9d5
ks_multi_himulti = ApproximateTwoSampleKSTest(df_pl_himulti.tduravg./df_pl_himulti.tdur_circ, df_pl_multi.tduravg./df_pl_multi.tdur_circ )

# ╔═╡ 4741c041-f315-4e45-b4cc-638ce16bbbcc
ad_multi_himulti = KSampleADTest(df_pl_himulti.tduravg./df_pl_himulti.tdur_circ, df_pl_multi.tduravg./df_pl_multi.tdur_circ )

# ╔═╡ 92c2f1dd-e4b6-4896-9ebe-1ffe3ee5f49c
(;p_KS_1m=pvalue(ks_single_multi), p_AD_1m = pvalue(ad_single_multi),
	p_KS_12=pvalue(ks_single_double),p_AD_12 = pvalue(ad_single_double), 
	p_KS_23=pvalue(ks_multi_himulti),p_AD_23 = pvalue(ad_multi_himulti))

# ╔═╡ adb97bb7-5a53-4583-a5fa-b4e8d521f6fb
ad_multi_3way = KSampleADTest(df_pl_himulti.tduravg./df_pl_himulti.tdur_circ, df_pl_double.tduravg./df_pl_double.tdur_circ, df_pl_single.tduravg./df_pl_single.tdur_circ )

# ╔═╡ a6ec4534-fada-4210-a780-2e124c4b4595
md"k-sample Anderson-Darling test:  p-value =  $(pvalue(ad_multi_3way))"

# ╔═╡ fe9ab16c-0eef-4061-b31f-281e74647174
begin
	median_teff = median(map(x->collect(x.teff)[1],eachrow(df_systems)))
	df_pl_loteff = df_pl |> @filter(_.teff<=median_teff) |> DataFrame;
	df_pl_hiteff = df_pl |> @filter(_.teff> median_teff) |> DataFrame;
	teff_kraft_break = 6200
	df_pl_hotterkraft = df_pl_no_teffcut |> @filter(_.teff>=teff_kraft_break) |> DataFrame;
	df_pl_coolerkraft = df_pl_no_teffcut |> @filter(_.teff< teff_kraft_break) |> DataFrame;
end;

# ╔═╡ f79b9f5e-1651-429b-8f6b-d6581a77707f
if make_plots
	plt_teff = plot(xlabel="Normalized Transit Duration", ylabel="Cumulative", legend=:topleft)
	#plot!(plt_teff,(x_plt),ntdur_all_cum.(x_plt), label="All",color=:green) 
		
	ntdur_loteff_cum = ecdf(df_pl_loteff.tduravg./df_pl_loteff.tdur_circ)
	plot!(plt_teff,(x_plt), ntdur_loteff_cum.(x_plt), label="$(round(Int64,min_teff)) - $(round(Int64,median_teff))K", color=:red)
	ntdur_hiteff_cum = ecdf(df_pl_hiteff.tduravg./df_pl_hiteff.tdur_circ)
	plot!(plt_teff,(x_plt), ntdur_hiteff_cum.(x_plt), label="$(round(Int64,median_teff)) - $(round(Int64,max_teff)) K", color=:blue)
	title!("Divided by median Teff ($(round(Int64,median_teff))K)")
	xlims!(0.0,2.0)
end

# ╔═╡ 42a5cd4e-130a-4f65-9698-f5bb8fc2b97d
if make_kraft_plot
	plt_kraft = plot(xlabel="Normalized Transit Duration", ylabel="Cumulative", legend=:topleft)
	#plot!(plt_kraft,(x_plt),ntdur_all_cum.(x_plt), label="All") 
		
	ntdur_coolerkraft_cum = ecdf(df_pl_coolerkraft.tduravg./df_pl_coolerkraft.tdur_circ)
	plot!(plt_kraft,(x_plt), ntdur_coolerkraft_cum.(x_plt), label="T_eff < 6200K", color=:red)
	ntdur_hotterkraft_cum = ecdf(df_pl_hotterkraft.tduravg./df_pl_hotterkraft.tdur_circ)
	
	plot!(plt_kraft,(x_plt), ntdur_hotterkraft_cum.(x_plt), label="T_eff > 6200K", color=:blue)
	title!(plt_kraft,"Divided by Kraft Break")
	xlims!(plt_kraft,0,2.0)
end

# ╔═╡ 89a0f47a-4b7f-43f6-b1fb-51db38c771b2
ad_teff = KSampleADTest(df_pl_loteff.tduravg./df_pl_loteff.tdur_circ, df_pl_hiteff.tduravg./df_pl_hiteff.tdur_circ)

# ╔═╡ 5b2daca0-997b-448c-98b4-690e9df6253b
ks_teff = ApproximateTwoSampleKSTest(df_pl_loteff.tduravg./df_pl_loteff.tdur_circ, df_pl_hiteff.tduravg./df_pl_hiteff.tdur_circ )

# ╔═╡ 0557efb5-6120-4c48-a9b4-3da2854d2f98
(;p_KS=pvalue(ks_teff),p_AD = pvalue(ad_teff))

# ╔═╡ d6b8477f-02a4-4c0b-ab80-4a1b0f2757c7
ks_kraft = ApproximateTwoSampleKSTest(df_pl_coolerkraft.tduravg./df_pl_coolerkraft.tdur_circ, df_pl_hotterkraft.tduravg./df_pl_hotterkraft.tdur_circ)

# ╔═╡ 1440edf5-38f7-4c58-96d8-8a6cefaed9c6
ad_kraft = KSampleADTest(df_pl_coolerkraft.tduravg./df_pl_coolerkraft.tdur_circ, df_pl_hotterkraft.tduravg./df_pl_hotterkraft.tdur_circ)

# ╔═╡ b3cb8a30-8165-4725-85a6-1026303d8932
(;p_KS=pvalue(ks_kraft),p_AD = pvalue(ad_kraft))

# ╔═╡ 7c75fdc7-2cbe-4f83-8d89-25cb3d328d11
begin
	teff_lifetime = 5400
	df_pl_hotter_lifetime = df_pl |> @filter(_.teff>=teff_lifetime) |> DataFrame;
	df_pl_cooler_lifetime = df_pl |> @filter(_.teff< teff_lifetime) |> DataFrame;
end;

# ╔═╡ 628d4417-a7f4-4958-9c46-c2ed93fd51cb
begin
	plt_lifetime = plot(xlabel="τ", ylabel="Cumulative", legend=:topleft)
	plot!(plt_teff,(x_plt),ntdur_all_cum.(x_plt), label="All") 
		
	ntdur_cooler_lifetime_cum = ecdf(df_pl_cooler_lifetime.tduravg./df_pl_cooler_lifetime.tdur_circ)
	plot!(plt_lifetime,(x_plt), ntdur_cooler_lifetime_cum.(x_plt), label="$(round(Int64,min_teff)) - $teff_lifetime K", color=:red)
	ntdur_hotter_lifetime_cum = ecdf(df_pl_hotter_lifetime.tduravg./df_pl_hotter_lifetime.tdur_circ)
	plot!(plt_lifetime,(x_plt), ntdur_hotter_lifetime_cum.(x_plt), label="$teff_lifetime - $(round(Int64,max_teff)) K", color=:blue)
	xlims!(plt_lifetime,0,2)
	title!(plt_lifetime,"Divided by $(teff_lifetime)K (Moorhead+ 2011)")
end

# ╔═╡ 6da264a6-d6b1-4bfd-85ce-fd80ff69ed09
ks_lifetime = ApproximateTwoSampleKSTest(df_pl_hotter_lifetime.tduravg./df_pl_hotter_lifetime.tdur_circ, df_pl_cooler_lifetime.tduravg./df_pl_cooler_lifetime.tdur_circ )

# ╔═╡ 533a81b8-51fa-4ea6-9dbb-ebfd5efe99c3
ad_lifetime = KSampleADTest(df_pl_cooler_lifetime.tduravg./df_pl_cooler_lifetime.tdur_circ, df_pl_hotter_lifetime.tduravg./df_pl_hotter_lifetime.tdur_circ)

# ╔═╡ b7b189f5-1ccc-4536-9022-af317dd35ca7
(;p_KS=pvalue(ks_lifetime),p_AD = pvalue(ad_lifetime))

# ╔═╡ 6deb2237-5e6d-4dac-bcca-c5193f5a4c61
begin
	median_rpl = median(df_pl.radius)
	df_pl_small_rp = df_pl |> @filter(_.radius<=median_rpl) |> DataFrame;
	df_pl_large_rp = df_pl |> @filter(_.radius> median_rpl) |> DataFrame;
end;

# ╔═╡ 69d8355b-ad3f-4ecb-9710-17b5846d6d1b
if make_plots
	plt_rp = plot(xlabel="τ", ylabel="Cumulative", legend=:topleft)
	#plot!(plt_rp,(x_plt),ntdur_all_cum.(x_plt), label="All") 
		
	ntdur_smallrp_cum = ecdf(df_pl_small_rp.tduravg./df_pl_small_rp.tdur_circ)
	plot!(plt_rp,(x_plt), ntdur_smallrp_cum.(x_plt), label="Small R_p")
	ntdur_largerp_cum = ecdf(df_pl_large_rp.tduravg./df_pl_large_rp.tdur_circ)
	plot!(plt_rp,(x_plt), ntdur_largerp_cum.(x_plt), label="Large R_p")
	title!(plt_rp,"Grouped by Planet Size (median = $(round(median_rpl,digits=2))R_⊕)")
	xlims!(plt_rp,0.0,2.0)
end

# ╔═╡ c02c95b0-454d-435a-8331-ebdb50d17414
ks_Rp = ApproximateTwoSampleKSTest(df_pl_large_rp.tduravg./df_pl_large_rp.tdur_circ, df_pl_small_rp.tduravg./df_pl_small_rp.tdur_circ )

# ╔═╡ 10e7a5ea-2980-4122-9ef1-385467ce08bb
ad_Rp = KSampleADTest(df_pl_large_rp.tduravg./df_pl_large_rp.tdur_circ, df_pl_small_rp.tduravg./df_pl_small_rp.tdur_circ )

# ╔═╡ 924fffcd-b1fa-455a-95d6-dbea0ebaf7d4
(;p_KS=pvalue(ks_Rp),p_AD = pvalue(ad_Rp))

# ╔═╡ f8a18459-2375-41e2-bc03-43662b42f6a6
begin
	median_P = median(df_pl.TTVperiod)
	split_P = median_P # 10
	df_pl_small_P = df_pl |> @filter(_.TTVperiod<=split_P) |> DataFrame;
	df_pl_large_P = df_pl |> @filter(_.TTVperiod> split_P) |> DataFrame;
end;

# ╔═╡ d25b4d73-fe48-4958-937b-8d06bb71852b
if make_plots
	plt_P = plot(xlabel="τ", ylabel="Cumulative", legend=:topleft)
	plot!(plt_P,(x_plt),ntdur_all_cum.(x_plt), label="All") 
	
	ntdur_smallP_cum = ecdf(df_pl_small_P.tduravg./df_pl_small_P.tdur_circ)
	plot!(plt_P,(x_plt), ntdur_smallP_cum.(x_plt), label="P<$(round(split_P*10)/10) d")
	ntdur_largeP_cum = ecdf(df_pl_large_P.tduravg./df_pl_large_P.tdur_circ)
	plot!(plt_P,(x_plt), ntdur_largeP_cum.(x_plt), label="P>$(round(split_P*10)/10) d")
	xlims!(plt_P,0.0,2.0)
	#title!(plt_P,@sprintf("Grouped by Period (median=%2.2f)",median_P))
end

# ╔═╡ 37fdce37-1188-4f67-bf18-8041b73e023c
ks_P = ApproximateTwoSampleKSTest(df_pl_large_P.tduravg./df_pl_large_P.tdur_circ, df_pl_small_P.tduravg./df_pl_small_P.tdur_circ )

# ╔═╡ c3555ff6-f87f-4cc7-ba11-519496969e33
pvalue(ks_P)

# ╔═╡ 7900c30c-200f-4a57-82d6-3261f56017e3
ad_P = KSampleADTest(df_pl_large_P.tduravg./df_pl_large_P.tdur_circ, df_pl_small_P.tduravg./df_pl_small_P.tdur_circ )

# ╔═╡ 93364c3e-3538-415e-99e5-5c018753970d
(;p_KS=pvalue(ks_P),p_AD = pvalue(ad_P))

# ╔═╡ f08261cf-6a76-47f9-b858-08ef68e1cce6
md"## Tinkering"

# ╔═╡ ea903194-4ad8-462a-8668-ac464eb14f3c
names(df_pl)

# ╔═╡ bf40c596-cd67-4247-bed4-11d53a843c22
median(df_pl.rhostar)

# ╔═╡ 9cb5e6d5-0853-45ab-b1a4-b588893fd143
let 
	plt = histogram(df_pl[!,"rhostar"], nbins = 100, label=:none)
	xlabel!(plt,"ρ (g/cm³)")
	ylabel!(plt,"Count")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.9.10"
DataFrames = "~1.2.2"
Distributions = "~0.25.28"
HypothesisTests = "~0.10.6"
LaTeXStrings = "~1.3.0"
Plots = "~1.23.6"
PlutoTest = "~0.2.0"
PlutoUI = "~0.7.19"
Query = "~1.0.0"
StatsBase = "~0.33.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "74147e877531d7c172f70b492995bc2b5ca3a843"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.10"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "794daf62dce7df839b8ed446fc59c68db4b5182f"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.3.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "cab6fd4d6a0fca4d7f1dcdc2a130884e6ae242c9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.28"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FilePathsBase]]
deps = ["Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "5440c1d26aa29ca9ea848559216e5ee5f16a8627"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.14"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Random", "Rmath", "Roots", "Statistics", "StatsBase"]
git-tree-sha1 = "dc9bb7abfa265e0cf030635315184a476a2dd5f3"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.10.6"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "19cb49649f8c41de7fea32d089d37de917b553da"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.0.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c8b8775b2f242c80ea85c83714c64ecfa3c53355"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "0d185e8c33401084cab546a756b387b15f76720c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.23.6"

[[PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "92b8ae1eee37c1b8f70d3a8fb6c3f2d81809a1c5"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.0"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[Query]]
deps = ["DataValues", "IterableTables", "MacroTools", "QueryOperators", "Statistics"]
git-tree-sha1 = "a66aa7ca6f5c29f0e303ccef5c8bd55067df9bbe"
uuid = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
version = "1.0.0"

[[QueryOperators]]
deps = ["DataStructures", "DataValues", "IteratorInterfaceExtensions", "TableShowUtils"]
git-tree-sha1 = "911c64c204e7ecabfd1872eb93c49b4e7c701f02"
uuid = "2aef5ad7-51ca-5a8f-8e88-e75cf067b44b"
version = "0.9.3"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "4c40dc61b51054bdb93536400420d73fdca6865e"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.7"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "def0718ddbabeb5476e51e5a43609bee889f285d"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.0"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "385ab64e64e79f0cd7cfcf897169b91ebbb2d6c8"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.13"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─88d46a9d-adda-4750-87e8-58c1b50176cd
# ╟─cb9f36fe-8bae-4e54-ae9c-d2fdc466c425
# ╟─69f7e8be-3541-4fa7-b74b-7ebfc4a0fefc
# ╟─92ef31db-67ae-4d4c-8c91-69734ee64d7a
# ╟─eadfde5c-9124-4c9f-87b6-a5c79d50d772
# ╟─dc9d19bc-cdcb-4a19-bbda-71203b0832a3
# ╟─f8fa6349-04be-4f6c-9465-85208dfc198d
# ╟─05572fad-3d61-4a99-8511-6915b0cda15a
# ╟─3d97e21f-9785-46d8-aa61-87a69fda8b19
# ╟─08aeb059-fe92-455d-a726-577329d68be7
# ╟─d80424df-bcfa-4bf1-805e-954fd4c9c25c
# ╟─77de3d40-428f-463e-88a7-3eb209b52747
# ╟─4c77a70d-c003-48eb-a683-359446950878
# ╟─1e682cd4-2ac0-43d2-98df-a054520fd5c5
# ╟─33d289e5-a205-4a82-a0e6-0b09d50a1f1f
# ╟─8dc05082-c2bd-4120-b7f0-499287953056
# ╟─4276feec-2e6f-4a89-89be-7a152a3833dc
# ╟─dac4d650-c8b7-4a24-8420-a0ebc56dccca
# ╟─da90b894-32da-4187-b1a1-6548bf5b5b12
# ╟─561c5fee-e103-4de8-ab47-e11f32a2a9d0
# ╟─56effba4-4537-4f1b-8225-9683668ca37a
# ╟─b5d1da7a-c483-42f7-9b64-e31bedf768ba
# ╟─d3b2e526-8838-460c-ad02-a6bde5c5ef4f
# ╠═7fe5f5af-2085-49ab-9a6c-990548ccb5f3
# ╟─5af5ef4f-a8a4-48b8-8111-c30765cf0ac0
# ╟─ab590e37-4c5b-4609-876a-c09ee92985c8
# ╟─a6074ab0-eacd-427b-971a-20ee9af02abb
# ╟─a6ec4534-fada-4210-a780-2e124c4b4595
# ╟─92c2f1dd-e4b6-4896-9ebe-1ffe3ee5f49c
# ╟─9f930a87-597d-4997-9819-6a41acdd53bb
# ╟─69d8355b-ad3f-4ecb-9710-17b5846d6d1b
# ╟─924fffcd-b1fa-455a-95d6-dbea0ebaf7d4
# ╟─95dfd011-b505-4c94-93f8-9968e5eb6cab
# ╟─616bb0a9-f6da-458a-bc25-ad3dbcd0d77b
# ╟─d25b4d73-fe48-4958-937b-8d06bb71852b
# ╟─93364c3e-3538-415e-99e5-5c018753970d
# ╠═4cf1d6c9-ff78-4056-be09-2279ccb3b54c
# ╟─3dfa1096-7625-409c-a3a3-24f3907d4ac5
# ╟─f8cc2102-eeb2-4b4c-a8d0-ec5c206a7f42
# ╟─df6285f5-9e8a-4dee-bbd3-64f1915ed5ba
# ╟─c859043d-d7ba-47cd-9035-319cdc200720
# ╟─f79b9f5e-1651-429b-8f6b-d6581a77707f
# ╟─0557efb5-6120-4c48-a9b4-3da2854d2f98
# ╟─42a5cd4e-130a-4f65-9698-f5bb8fc2b97d
# ╟─b3cb8a30-8165-4725-85a6-1026303d8932
# ╟─628d4417-a7f4-4958-9c46-c2ed93fd51cb
# ╟─b7b189f5-1ccc-4536-9022-af317dd35ca7
# ╟─4ee2936e-d41f-41f1-b307-60d606629b5b
# ╟─98a95e8f-bc1d-4a2a-9665-0e661673b025
# ╠═89a0f47a-4b7f-43f6-b1fb-51db38c771b2
# ╠═5b2daca0-997b-448c-98b4-690e9df6253b
# ╠═6da264a6-d6b1-4bfd-85ce-fd80ff69ed09
# ╠═533a81b8-51fa-4ea6-9dbb-ebfd5efe99c3
# ╠═d6b8477f-02a4-4c0b-ab80-4a1b0f2757c7
# ╠═1440edf5-38f7-4c58-96d8-8a6cefaed9c6
# ╠═5cf370f9-875f-4298-bcb2-7481473b2a41
# ╠═860b0cdc-67d3-47ce-856f-10f172a97345
# ╠═ea72041d-6bdf-4263-9a78-884edaf2c0af
# ╟─ff82c99c-b489-43a8-abc4-ceaf694aabf8
# ╠═18dbba04-1be6-4a59-9d11-c2c95982897b
# ╠═15f7a4d8-d827-441a-b599-902b57f19821
# ╠═9a076cb1-f5fa-45fd-a624-071c2677691d
# ╠═2cd616cd-778e-490a-9826-343cf3b42c9d
# ╠═ad551ae3-e000-4259-98ea-4cad99f6b9d5
# ╠═4741c041-f315-4e45-b4cc-638ce16bbbcc
# ╠═adb97bb7-5a53-4583-a5fa-b4e8d521f6fb
# ╟─639dfbeb-60e5-4f59-a3dc-81a0884f0ac4
# ╠═c02c95b0-454d-435a-8331-ebdb50d17414
# ╠═10e7a5ea-2980-4122-9ef1-385467ce08bb
# ╠═bda80d86-3430-45ef-9006-d6b77b1758d8
# ╟─ae6f7834-f939-4a09-9bd6-bf105c6b05c5
# ╠═37fdce37-1188-4f67-bf18-8041b73e023c
# ╠═7900c30c-200f-4a57-82d6-3261f56017e3
# ╠═c3555ff6-f87f-4cc7-ba11-519496969e33
# ╠═3013ec6b-3139-4fdd-aa69-fe4b45dfd9f8
# ╟─d498d791-9763-47e5-b5cc-3a41d930dbe7
# ╠═893a58a5-f88b-4a59-bf21-c66b7e9404d4
# ╠═46749afb-7cdb-4956-b835-8623b7cd856a
# ╠═bc6d7cfc-b8da-43ac-84a5-73e323323701
# ╠═cd152b2a-1690-42d5-b76a-cd71b636aa01
# ╟─7eb040e9-300b-41d2-9171-9f48c9b98758
# ╠═55fcfaca-03b9-4d82-822c-5ddf22bc350e
# ╠═61c1412d-7773-4a26-ba62-c3e91ff5d908
# ╠═859b2e3f-7c97-4d3f-b3a0-cb19f5da6a8a
# ╠═c7fc5e4b-c1bc-4e65-a126-34b9dacce8c3
# ╠═33a32461-9f78-48e0-8e33-9cf767c18d9b
# ╠═49c7d182-fd27-4d51-a7d2-c2153beda6af
# ╠═980de46b-4839-4fd7-b4fe-4ebe1b19bf43
# ╟─2704c820-d307-4acc-869a-e8670232ec38
# ╠═35d1e842-496b-11ec-2e35-17c2d536f22f
# ╠═cc668221-1b44-4fcb-83bf-d7f851c8bcbd
# ╟─585c1c7a-56bd-49b7-b075-b2232072c893
# ╠═3a0a2bd4-32f4-47d7-a34d-b22880c91395
# ╠═98e17af2-f0b4-4a2a-a9ba-2040da85af88
# ╠═e4581234-6052-4ec4-b6d1-94d3c3c10a28
# ╠═cc173a8d-09a4-4ff0-a695-6eba16ade3b1
# ╠═ac1dcf04-19e5-4221-891f-99747b9cf5f4
# ╠═5eccb812-d7bb-4f9f-a7a0-6892d8eb46ac
# ╠═1a7dc30e-a793-435c-9c79-c96ad455bb11
# ╠═74fd8933-74f4-4364-a2f8-1465740987b9
# ╟─591cebe7-9aa4-4617-9aee-c10c038a7002
# ╠═0ac069bb-3e86-4648-964e-69c640538187
# ╠═ed49e14f-4b9c-42d2-8e55-2b3091fe4c14
# ╠═9179cdfd-1f7a-476f-a535-bae6a2fc192f
# ╠═b3447847-e7f2-4844-be44-01669a0f0e74
# ╠═8f46b232-4978-4ea5-bade-d8446e37b8ff
# ╠═2f44dc3d-00f9-4fb9-b4af-e15eead627c2
# ╠═52ca81e0-fa57-4957-a91d-4b6919c3366e
# ╠═0c3d9129-c9be-4dc4-a1a0-1b8ec39ba16c
# ╠═5d9c5e6e-d42e-4bf0-99a0-fc89612c15d8
# ╟─9dbf4174-17c6-421f-ba9e-953221fcd58f
# ╠═b97166df-bf9a-4208-ac69-6475e6e1cdfa
# ╠═fe9ab16c-0eef-4061-b31f-281e74647174
# ╠═7c75fdc7-2cbe-4f83-8d89-25cb3d328d11
# ╠═6deb2237-5e6d-4dac-bcca-c5193f5a4c61
# ╠═f8a18459-2375-41e2-bc03-43662b42f6a6
# ╟─f08261cf-6a76-47f9-b858-08ef68e1cce6
# ╠═ea903194-4ad8-462a-8668-ac464eb14f3c
# ╠═bf40c596-cd67-4247-bed4-11d53a843c22
# ╠═9cb5e6d5-0853-45ab-b1a4-b588893fd143
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
