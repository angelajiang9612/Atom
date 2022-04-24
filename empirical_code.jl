
using Pkg
using DataFrames
using CSV
using Statistics
using DataStructures
using Distributions
using DelimitedFiles
using Random
cd("/Users/bubbles/Desktop/751 Labor/751 Chao/Empirical Project")
data=readdlm("data_age4554.txt")
df = DataFrames.DataFrame(data, :auto)
df = df[2:30121,:]
df=rename!(df,[:id,:age,:lfp,:x,:wage,:educ,:lfp0,:hinc])

K=45 ##the maximum first period experience level, which is 24, +20, K is length so add 1
k_grid=Int.(collect(range(0, stop = 44, length = K)))
K_max=K-1
A_final=64 ##A for terminating period
A_init=45
T_data=10 ##data periods
δ=0.95

function working_by_age(df_data)
    participation_by_age=zeros(10)
    for t=45:54
        df_age=filter(row -> row.age==t, df_data)
        participation_by_age[t-44]=mean(df_age.lfp)
    end
    return participation_by_age
end


function wages_by_age(df_data_working)
    wages_by_age=zeros(10)
    for t=45:54
        df_age=filter(row -> row.age==t, df_data_working)
        wages_by_age[t-44]=mean(df_age.wage)
    end
    return wages_by_age
end


#5 labor participation by lagged labor participation

function L_estimate(df,l) #this function calculates the destination firm types for those leaving from type l
    counter_0=0  ##counters
    counter_1=0
    for i=1:N
        i_rows=Matrix(df[df.id .== i, :])
        i_work=i_rows[:,3] ##
        indexes=findall(lfp->lfp==l, i_work)
        for t in indexes ##only look at those with labor=l
            if t < 10 #if not the last period
                new_lfp=Integer(i_work[t+1])
                if new_lfp == 0
                    counter_0=counter_0+1
                elseif new_lfp ==1
                    counter_1=counter_1+1
                end
            end
        end
    end
    return counter_0, counter_1
end

transite_00, transite_01 =L_estimate(df_data,0) ##people in this data seem to always work or always not work
transite_10, transite_11 =L_estimate(df_data,1)


Transition_Matrix=[transite_00/(transite_00 +transite_01) transite_01/(transite_00 +transite_01); transite_10/(transite_10+transite_11) transite_11/(transite_10+transite_11)]



function cutoffs_analytical() ##for akk the index k, need to add 1
    E_MAX=zeros(N,K,T_data)
    xi_star_matrix=zeros(N,K,T_data)  ##place to store cutoffs
    for i=1:N
        println("this is", i, "person")
        for t=T_data:-1:1
            row=Matrix(filter(row -> row.id ==i && row.age==44+t, df))
            y=row[8]
            s=row[6]
            if t==T_data #for last period
                for k=0:K_max #last period possible experience 0 to k_max, k is experience
                    xi_star=log(-α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s) - (β_1 + β_2*k + β_3*k^2 + β_4*s) - log(1 +α_2) #xi_star us calculated for any state variable k experience in the last year
                    X_A= exp(β_1 + β_2*k + β_3*k^2 + β_4*s)
                    xi_star_matrix[i,k+1,t]=xi_star #to save need to add one to k to get index starting from. Saving e*(k) for all possible k
                    E_MAX[i,k+1,t]=y*(1-cdf.(Normal(),xi_star/σ_ξ)) + X_A*exp(0.5*σ_ξ^2)*(1-cdf.(Normal(),(xi_star-σ_ξ^2)/σ_ξ)) + ((1+α_2)*y + α_1)*cdf.(Normal(),xi_star/σ_ξ) #E_T, for use in T-1
                end
            elseif t<=T_data-1
                for k=0:K_max-(T_data-t)
                    xi_star=log((-α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s) + δ*(E_MAX[i,k+1,t+1]-E_MAX[i,k+1+1,t+1])) - (β_1 + β_2*k + β_3*k^2 + β_4*s) - log(1+α_2) #T-1 period cutoff only depends on expected for T period
                    X_A= exp(β_1 + β_2*k + β_3*k^2 + β_4*s)
                    E_MAX[i,k+1,t]=(y + δ*E_MAX[i,k+1+1,t+1])*(1-cdf.(Normal(),xi_star/σ_ξ)) + X_A*exp(0.5*σ_ξ^2)*(1-cdf.(Normal(),(xi_star-σ_ξ^2)/σ_ξ)) +
                    ((1+α_2)*y + α_1 + E_MAX[i,k+1,t+1])*cdf.(Normal(),xi_star/σ_ξ) #this E_mAx now has uncertainty in both work and not work world.
                    xi_star_matrix[i,k+1,t]=xi_star
                end
            end
        end
    end
    return xi_star_matrix
end


function getendo_new(xi_star_matrix, xi_matrix) ##this generates endogenous x and k
    k_vector=zeros(N,T_data,M)
    x=zeros(N,T_data,M)
    for m=1:M
        for i=1:N
            println("this is", i, "person")
            row=Matrix(filter(row -> row.id ==i && row.age==45, df))
            for t=1:T_data
                xi=xi_matrix[i,t,m]
                if t==1
                    k_vector[i,t,m]=row[4]
                    k= row[4]
                    xi_cutoff=xi_star_matrix[i,k+1,t]
                    if xi >= xi_cutoff
                        k_vector[i,t+1,m]=k+1 #update experience next period
                        x[i,t,m]=1 #t
                    elseif xi < xi_cutoff
                        k_vector[i,t+1,m]=k
                        x[i,t,m]=0
                    end
                elseif 1<t<T_data
                    k= floor(Int,k_vector[i,t,m])
                    xi_cutoff=xi_star_matrix[i,k+1,t]
                    if xi >= xi_cutoff
                        k_vector[i,t+1,m]=k+1 #update experience next period
                        x[i,t,m]=1 #this period participation is 1
                    elseif xi < xi_cutoff
                        k_vector[i,t+1,m]=k
                        x[i,t,m]=0
                    end
                else t==T_data #do not need to update next period.
                    k= floor(Int,k_vector[i,t,m])
                    xi_cutoff=xi_star_matrix[i,k+1,t]
                    if xi >= xi_cutoff
                        x[i,t,m]=1 #this period participation is 1
                    elseif xi < xi_cutoff
                        x[i,t,m]=0
                    end
                end
            end
        end
    end
    return x, k_vector
end


function simulate_data(lfp_vector,k_vector,wages_vector,df_for_calib) ##create function to store all sim_data
    s_id=zeros(N*T_data)
    s_age=zeros(N*T_data)
    s_lfp=zeros(N*T_data)
    s_x=zeros(N*T_data)
    s_wages=zeros(N*T_data)
    s_educ=zeros(N*T_data)
    sim_data=zeros(N*T_data,6,M)

    for m=1:M
        s_id=df_for_calib.id
        s_age=df_for_calib.age
        s_lfp=vec(lfp_vector[:,:,m]') #need ' to go by row
        s_x=vec(k_vector[:,:,m]')
        s_wages=vec(wages_vector[:,:,m]')
        s_educ=df_for_calib.educ
        sim=DataFrames.DataFrame(hcat(s_id, s_age,s_lfp, s_x,s_wages,s_educ), :auto)
        sim_data[:,:,m].=rename!(sim,[:id,:age,:lfp,:x,:wage,:educ])
    end
    return sim_data
end


##Parameters 

N=1506
M=2 ##S is number of simulations
#set initial parameters
##sigma eta is randomly set (=rho), b is randomly set, the rest is using results from Wholpin paper
α_1=-1375
α_2=-0.047
α_3=-11.5
α_4=-95.6
b=0.1
β_1=-0.280
β_2=0.024
β_3=-0.0002
β_4=0.050
σ_ξ=0.194



function get_moments_difference()
    Random.seed!(1234);
    xi_matrix = reshape(rand(Normal(0, σ_ξ), N*T_data*M), N,T_data,M)
    xi_star_matrix= cutoffs_analytical()
    lfp_vector, k_vector = getendo_new(xi_star_matrix, xi_matrix)
    df_for_calib=filter(row -> row.age<=54 && row.id<=N, df)
    educ=reshape(df_for_calib.educ, N,T_data)
    wages_vector = exp.(β_1.+ β_2.*k_vector.+ β_3.*k_vector.^2 .+ xi_matrix.+ educ) ##need to check this, this generates wages
    sim_data_all=simulate_data(lfp_vector,k_vector,wages_vector,df_for_calib)

    working_by_educ_diff=zeros(5,M)
    working_by_age_diff=zeros(10,M)
    working_by_experience_diff=zeros(3,M)
    wages_by_educ_diff=zeros(5,M)
    wages_by_age_diff=zeros(10,M)
    transition_diff=zeros(2,2,M)

    
    df_for_calib=filter(row -> row.age<=54 && row.id<=N, df)  #data stuff that doesn't depend on simulation
    df_data_educ_under_12=filter(row -> row.educ<=11, df_for_calib)
    df_data_educ_12=filter(row -> row.educ==12, df_for_calib)
    df_data_educ_13_15=filter(row -> row.educ<=15,df_for_calib)
    df_data_educ_16_plus=filter(row -> row.educ>=16, df_for_calib)
    df_data_exp_10_or_less=filter(row -> row.x<=10, df_for_calib)
    df_data_exp_11_to_20=filter(row -> 11<=row.x<=20, df_for_calib)
    df_data_exp_21_plus=filter(row -> row.x >= 21,  df_for_calib)

    df_data_working=filter(row -> row.lfp==1, df_for_calib)
    df_data_working_educ_under_12=filter(row -> row.educ<=12, df_data_working)
    df_data_working_educ_12=filter(row -> 11<=row.educ==12, df_data_working)
    df_data_working_educ_13_15=filter(row -> 13 <= row.educ <= 15,  df_data_working)
    df_data_working_educ_16_plus=filter(row -> row.educ>=16,  df_data_working)

    transite_00_data, transite_01_data =L_estimate(df_for_calib,0) ##people in this data seem to always work or always not work
    transite_10_data, transite_11_data =L_estimate(df_for_calib,1)
    Transition_Matrix=[transite_00_data/(transite_00_data +transite_01_data) transite_01_data/(transite_00_data +transite_01_data); transite_10_data/(transite_10_data+transite_11_data) transite_11_data/(transite_10_data+transite_11_data)]


    for m=1:M
        sim_data=DataFrames.DataFrame(sim_data_all[:,:,m], :auto)
        sim_data=rename!(sim_data,[:id,:age,:lfp,:x,:wage,:educ])
        sim_data_educ_under_12=filter(row -> row.educ<=11, sim_data)
        sim_data_educ_12=filter(row -> row.educ==12, sim_data)
        sim_data_educ_13_15=filter(row -> row.educ<=15, sim_data)
        sim_data_educ_16_plus=filter(row -> row.educ>=16, sim_data)
        sim_data_working=filter(row -> row.lfp==1, sim_data)
        sim_data_exp_10_or_less=filter(row -> row.x<=10, sim_data)
        sim_data_exp_11_to_20=filter(row -> 11<=row.x<=20, sim_data)
        sim_data_exp_21_plus=filter(row -> row.x >= 21, sim_data)
        sim_data_working_educ_under_12=filter(row -> row.educ<=12, sim_data_working)
        sim_data_working_educ_12=filter(row -> 11<=row.educ==12, sim_data_working)
        sim_data_working_educ_13_15=filter(row -> 13 <= row.educ <= 15,  sim_data_working)
        sim_data_working_educ_16_plus=filter(row -> row.educ>=16,  sim_data_working)

        working_by_educ_diff[1,m]=mean(df_for_calib.lfp)-mean(sim_data.lfp)
        working_by_educ_diff[2,m]=mean(df_data_educ_under_12.lfp)-mean(sim_data_educ_under_12.lfp)
        working_by_educ_diff[3,m]=mean(df_data_educ_12.lfp)-mean(sim_data_educ_12.lfp)
        working_by_educ_diff[4,m]=mean(df_data_educ_13_15.lfp)- mean(sim_data_educ_13_15.lfp)
        working_by_educ_diff[5,m]=mean(df_data_working_educ_16_plus.lfp)-mean(sim_data_educ_16_plus.lfp)

        participation_by_age=working_by_age(df_for_calib)
        sim_participation_by_age=working_by_age(sim_data)

        working_by_age_diff[:,m]=participation_by_age.-sim_participation_by_age

        working_by_experience_diff[:,m]=[mean(df_data_exp_10_or_less.lfp)-mean(sim_data_exp_10_or_less.lfp) mean(df_data_working_exp_11_to_20.lfp)- mean(sim_data_exp_11_to_20.lfp) mean(df_data_working_exp_21_plus.lfp)- mean(sim_data_exp_21_plus.lfp)
        ]


        wages_by_educ_diff[:,m]=[mean(df_data_working.wage)-mean(sim_data_working.wage) mean(df_data_working_educ_under_12.lfp)-mean(sim_data_working_educ_under_12.lfp) mean(df_data_working_educ_12.lfp)- mean(sim_data_working_educ_12.lfp) mean(df_data_working_educ_13_15.lfp)- mean(sim_data_working_educ_13_15.lfp) mean(df_data_working_educ_16_plus.lfp)-mean(sim_data_working_educ_16_plus.lfp)]
        wages_age=wages_by_age(df_data_working)
        sim_wages_age=wages_by_age(sim_data_working)
        wages_by_age_diff[:,m]=wages_age.-sim_wages_age

        transite_00_sim, transite_01_sim =L_estimate(sim_data,0) ##people in this data seem to always work or always not work
        transite_10_sim, transite_11_sim =L_estimate(sim_data,1)
        Transition_Matrix_sim=[transite_00_sim/(transite_00_sim +transite_01_sim) transite_01_sim/(transite_00_sim +transite_01_sim); transite_10_sim/(transite_10_sim+transite_11_sim) transite_11_sim/(transite_10_sim+transite_11_sim)]

        transition_diff[:,:,m]=Transition_Matrix.-Transition_Matrix_sim
    end
        working_by_educ_diff_mean=mean(working_by_educ_diff,dims=2)
        working_by_age_diff_mean=mean(working_by_age_diff,dims=2)
        working_by_experience_diff_mean=mean(working_by_experience_diff,dims=2)
        wages_by_educ_diff_mean=mean(wages_by_educ_diff,dims=2)
        wages_by_age_diff_mean=mean(wages_by_age_diff,dims=2)
        transition_diff_mean=mean(transition_diff,dims=3)

        return working_by_educ_diff_mean, working_by_age_diff_mean, working_by_experience_diff_mean, wages_by_educ_diff_mean, wages_by_age_diff_mean, transition_diff_mean
end



@elapsed m1,m2,m3,m4,m5,m6 =get_moments_difference()

