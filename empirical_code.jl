using DataFrames
using CSV
using Statistics
using DataStructures
using Distributions
using DelimitedFiles
using Random


#importing data

cd("/Users/bubbles/Desktop/751 Labor/751 Chao/Empirical Project")
data=readdlm("data_age4554.txt")
df = DataFrames.DataFrame(data, :auto)
df = df[2:30121,:]
df=rename!(df,[:id,:age,:lfp,:x,:wage,:educ,:lfp0,:hinc])


##Parameters that do not change

K=35 ##the maximum first period experience level, which is 24, +10 for now, K is length so add 1. Don't need to go any higher at this point
k_grid=Int.(collect(range(0, stop = 34, length = K)))
K_max=K-1
A_final=64 ##A for terminating period
A_init=45
T_data=10 ##data periods
δ=0.95


#The three function below calculate moments from real data/simulated data, they should be fine because we have checked our results

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


function cutoffs_analytical(xi_matrix) ####This functionmcalculates cutoffs
    E_MAX=zeros(N,K,T_data)
    xi_star_matrix=zeros(N,K,T_data)  ##place to store cutoffs. This function use real data rather than grid so don't need grid for the other state variables
    for i=1:N
        println("this is", i, "person")
        for t=T_data:-1:1 #from backwards
            row=Matrix(filter(row -> row.id ==i && row.age==44+t, df)) #select i's info in time t. e.g.if T_data is 10, ag =44+10=54, which is the last period in data
            y=row[8] #the last row is husband's income
            s=row[6] #the sixth row is education
            if t==T_data #for last period
                for k=0:K_max #last period possible experience 0 to k_max of 34, k is experience so start from 0
                    log_inside= -α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s #dealing with the domain error problem
                    if log_inside >= 0
                        xi_star=log(-α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s) - log(1 +α_2)- (β_1 + β_2*k + β_3*k^2 + β_4*s)  #xi_star us calculated for any state variable k experience in the last year
                        xi_star_matrix[i,k+1,t]=xi_star #to save need to add one to k to get index.
                    elseif log_inside <0
                        xi_star=-1000000000000
                        xi_star_matrix[i,k+1,t]=xi_star
                    end
                    X_A= exp(β_1 + β_2*k + β_3*k^2 + β_4*s)
                    E_MAX[i,k+1,t]=(α_1 + (1+α_2)*(y-b) + α_3*k + α_4*s)*(1-cdf.(Normal(),xi_star/σ_ξ)) + (1+α_2)*X_A*exp(0.5*σ_ξ^2)*(1-cdf.(Normal(),(xi_star-σ_ξ^2)/σ_ξ)) + y*cdf.(Normal(),xi_star/σ_ξ)
                end
            elseif t<=T_data-1 #if not last period
                for k=0:K_max-(T_data-t) #e.g. if this period is t=6, the max experience possible is 34-(10-6)=30
                    log_inside = -α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s + δ*(E_MAX[i,k+1,t+1]-E_MAX[i,k+2,t+1])
                    if log_inside >= 0 #dealing with the domain error problem
                        xi_star=log(-α_1 - α_2*y + b*(1+α_2) - α_3*k -α_4*s + δ*(E_MAX[i,k+1,t+1]-E_MAX[i,k+2,t+1])) - log(1+α_2)-(β_1 + β_2*k + β_3*k^2 + β_4*s)  #T-1 period cutoff only depends on expected for T period
                        xi_star_matrix[i,k+1,t]=xi_star
                    elseif log_inside <0
                        xi_star=-1000000000000
                        xi_star_matrix[i,k+1,t]=xi_star
                    end
                    X_A= exp(β_1 + β_2*k + β_3*k^2 + β_4*s)
                    E_MAX[i,k+1,t]=(α_1 + (1+α_2)*(y-b) + α_3*k + α_4*s + E_MAX[i,k+1+1,t+1])*(1-cdf.(Normal(),xi_star/σ_ξ)) + (1+α_2)*X_A*exp(0.5*σ_ξ^2)*(1-cdf.(Normal(),(xi_star-σ_ξ^2)/σ_ξ))
                    + (y + E_MAX[i,k+1,t+1])*cdf.(Normal(),xi_star/σ_ξ)
                end
            end
        end
    end
    return xi_star_matrix
end


function getendo_new(xi_star_matrix, xi_matrix) ##this generates endogenous x and k for each simulation using initial experience
    k_vector=zeros(N,T_data,M) #experience of each person at each period in each simulation
    x=zeros(N,T_data,M) #action of each person at each period in each simulation
    for m=1:M
        for i=1:N
            println("this is", i, "person")
            row=Matrix(filter(row -> row.id ==i && row.age==45, df)) #only need the starting experience
            for t=1:T_data
                xi=xi_matrix[i,t,m] #get shock from shock matrix
                if t==1
                    k=row[4]
                    k_vector[i,t,m]=k #experience is 4th column, need to make sure this is an integer
                    xi_cutoff=xi_star_matrix[i,k+1,t] #add 1 for matrix indice
                    if xi >= xi_cutoff
                        x[i,t,m]=1 #t
                        k_vector[i,t+1,m]=k+1 #update experience next period
                    elseif xi < xi_cutoff
                        k_vector[i,t+1,m]=k
                        x[i,t,m]=0
                    end
                elseif 1<t<T_data
                    k=floor(Int,k_vector[i,t,m]) ##I don't know why I need to convert as interger again but if I don't do it I get error.
                    xi_cutoff=xi_star_matrix[i,k+1,t]
                    if xi >= xi_cutoff
                        k_vector[i,t+1,m]=k+1 #update experience next period
                        x[i,t,m]=1 #this period participation is 1
                    elseif xi < xi_cutoff
                        k_vector[i,t+1,m]=k
                        x[i,t,m]=0
                    end
                else t==T_data #do not need to update next period.
                    k=floor(Int,k_vector[i,t,m])
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



function simulate_data(lfp_vector,k_vector,wages_vector,df_for_calib) ##create function to store all simulated data in a datafram similar to the real data
    s_id=zeros(N*T_data) #create vectors store information, the number of rows is N times time period
    s_age=zeros(N*T_data)
    s_lfp=zeros(N*T_data)
    s_x=zeros(N*T_data)
    s_wages=zeros(N*T_data)
    s_educ=zeros(N*T_data)
    sim_data=zeros(N*T_data,6,M) ##to store the six things we care about, no initial labor or husband's wage
    for m=1:M
        s_id=df_for_calib.id #s_id is just id information from the df
        s_age=df_for_calib.age
        s_lfp=vec(lfp_vector[:,:,m]') #need ' to go by id first go through all t
        s_x=vec(k_vector[:,:,m]')
        s_wages=vec(wages_vector[:,:,m]')
        s_educ=df_for_calib.educ
        sim=DataFrames.DataFrame(hcat(s_id, s_age,s_lfp, s_x,s_wages,s_educ), :auto) #convert to dataframe format
        sim_data[:,:,m].=rename!(sim,[:id,:age,:lfp,:x,:wage,:educ]) #rename data
    end
    return sim_data
end


################Simulation part#############

N=1506
M=1 ##S is number of simulations
#set initial parameters
##sigma eta is randomly set (=rho), b is randomly set, the rest is using results from Wholpin paper
α_1=-14500
α_2=-0.25 #this cannot be less than -1 because we will get domain error
α_3=-100
α_4=-200
b=1000

β_1=9.93
β_2=0.01 #this term should also be small because wages doesnt change much with age, which means probably won't change much with experience 
β_3=-0.0002
β_4=0.035 #this term should be small because wage doesn't change much as education changes 
σ_ξ=0.194

#higher b1 and lower b4 might be needed for education fit 

@elapsed m1,m2,m3,m4,m5,m6 =get_moments_difference()
@show m1,m2,m3,m4,m5,m6


xi_star_matrix=cutoffs_analytical()

function get_moments_difference() #this compute moment difference

    df_for_calib=filter(row -> row.age<=54 && row.id<=N, df)  #real data stuff that doesn't depend on simulation but depend on N
    df_data_educ_under_12=filter(row -> row.educ<=11, df_for_calib)
    df_data_educ_12=filter(row -> row.educ==12, df_for_calib)
    df_data_educ_13_15=filter(row -> 13<=row.educ<=15,df_for_calib)
    df_data_educ_16_plus=filter(row -> row.educ>=16, df_for_calib)
    df_data_exp_10_or_less=filter(row -> row.x<=10, df_for_calib)
    df_data_exp_11_to_20=filter(row -> 11<=row.x<=20, df_for_calib)
    df_data_exp_21_plus=filter(row -> row.x >= 21,  df_for_calib)
    df_data_working=filter(row -> row.lfp==1, df_for_calib)
    df_data_working_educ_under_12=filter(row -> row.educ<=12, df_data_working)
    df_data_working_educ_12=filter(row -> row.educ==12, df_data_working)
    df_data_working_educ_13_15=filter(row -> 13 <= row.educ <= 15,  df_data_working)
    df_data_working_educ_16_plus=filter(row -> row.educ>=16,  df_data_working)
    transite_00_data, transite_01_data =L_estimate(df_for_calib,0) ##people in this data seem to always work or always not work
    transite_10_data, transite_11_data =L_estimate(df_for_calib,1)
    Transition_Matrix=[transite_00_data/(transite_00_data +transite_01_data) transite_01_data/(transite_00_data +transite_01_data); transite_10_data/(transite_10_data+transite_11_data) transite_11_data/(transite_10_data+transite_11_data)]
#question, do I have to write as function of xi_star etc all the time?
    Random.seed!(1234); #generating simulated results
    xi_matrix = reshape(rand(Normal(0, σ_ξ), N*T_data*M), N,T_data,M) #draw shocks
    xi_star_matrix= cutoffs_analytical(xi_matrix) #compute cutoffs
    lfp_vector, k_vector = getendo_new(xi_star_matrix, xi_matrix) #this calculates participation and experience
    df_for_calib=filter(row -> row.age<=54 && row.id<=N, df) #when I try this for small N I need to make sure I subset the real data by small N too
    educ_change=reshape(df_for_calib.educ,T_data,N) #need to get pick up education to calculate wages
    educ=educ_change' #do this to get in the right form
    wages_vector = exp.(β_1.+ β_2.*k_vector.+ β_3.*k_vector.^2 .+  β_4*educ.+xi_matrix)

    sim_data_all=simulate_data(lfp_vector,k_vector,wages_vector,df_for_calib) #convert the results for all the simulations into a datraframe, this is in N*T, 6, m dimension

    working_by_educ_diff=zeros(5,M) ##setting up matrix for storing the moment differences
    working_by_age_diff=zeros(10,M)
    working_by_experience_diff=zeros(3,M)
    wages_by_educ_diff=zeros(5,M)
    wages_by_age_diff=zeros(10,M)
    transition_diff=zeros(2,2,M)

    for m=1:M
        sim_data=DataFrames.DataFrame(sim_data_all[:,:,m], :auto) #take out one simulation's data
        sim_data=rename!(sim_data,[:id,:age,:lfp,:x,:wage,:educ])
        sim_data_educ_under_12=filter(row -> row.educ<=11, sim_data)
        sim_data_educ_12=filter(row -> row.educ==12, sim_data)
        sim_data_educ_13_15=filter(row -> 13<=row.educ<=15,sim_data)
        sim_data_educ_16_plus=filter(row -> row.educ>=16, sim_data)
        sim_data_exp_10_or_less=filter(row -> row.x<=10, sim_data)
        sim_data_exp_11_to_20=filter(row -> 11<=row.x<=20, sim_data)
        sim_data_exp_21_plus=filter(row -> row.x >= 21, sim_data)

        sim_data_working=filter(row -> row.lfp==1, sim_data)
        sim_data_working_educ_under_12=filter(row -> row.educ<=11, sim_data_working)
        sim_data_working_educ_12=filter(row -> row.educ==12, sim_data_working)
        sim_data_working_educ_13_15=filter(row -> 13 <= row.educ <= 15,  sim_data_working)
        sim_data_working_educ_16_plus=filter(row -> row.educ>=16,  sim_data_working)

        participation_by_age=working_by_age(df_for_calib)
        sim_participation_by_age=working_by_age(sim_data)

        working_by_educ_diff[1,m]=mean(df_for_calib.lfp)-mean(sim_data.lfp) #for each simulation compare results and store
        working_by_educ_diff[2,m]=mean(df_data_educ_under_12.lfp)-mean(sim_data_educ_under_12.lfp)
        working_by_educ_diff[3,m]=mean(df_data_educ_12.lfp)-mean(sim_data_educ_12.lfp)
        working_by_educ_diff[4,m]=mean(df_data_educ_13_15.lfp)- mean(sim_data_educ_13_15.lfp)
        working_by_educ_diff[5,m]=mean(df_data_educ_16_plus.lfp)-mean(sim_data_educ_16_plus.lfp)

        working_by_age_diff[:,m]=participation_by_age.-sim_participation_by_age

        working_by_experience_diff[:,m]=[mean(df_data_exp_10_or_less.lfp)-mean(sim_data_exp_10_or_less.lfp) mean(df_data_exp_11_to_20.lfp)-
        mean(sim_data_exp_11_to_20.lfp) mean(df_data_exp_21_plus.lfp)- mean(sim_data_exp_21_plus.lfp)]

        wages_by_educ_diff[1,m]=mean(df_data_working.wage)-mean(sim_data_working.wage)
        wages_by_educ_diff[2,m]=mean(df_data_working_educ_under_12.wage)-mean(sim_data_working_educ_under_12.wage)
        wages_by_educ_diff[3,m]=mean(df_data_working_educ_12.wage)- mean(sim_data_working_educ_12.wage)
        wages_by_educ_diff[4,m]= mean(df_data_working_educ_13_15.wage)- mean(sim_data_working_educ_13_15.wage)
        wages_by_educ_diff[5,m]=mean(df_data_working_educ_16_plus.wage)-mean(sim_data_working_educ_16_plus.wage)

        wages_age=wages_by_age(df_data_working)
        sim_wages_age=wages_by_age(sim_data_working)
        wages_by_age_diff[:,m]=wages_age.-sim_wages_age

        transite_00_sim, transite_01_sim =L_estimate(sim_data,0) ##people in this data seem to always work or always not work
        transite_10_sim, transite_11_sim =L_estimate(sim_data,1)
        Transition_Matrix_sim=[transite_00_sim/(transite_00_sim +transite_01_sim) transite_01_sim/(transite_00_sim +transite_01_sim);
        transite_10_sim/(transite_10_sim+transite_11_sim) transite_11_sim/(transite_10_sim+transite_11_sim)]

        transition_diff[:,:,m]=Transition_Matrix.-Transition_Matrix_sim
    end
        working_by_educ_diff_mean=mean(working_by_educ_diff,dims=2) #average over the second dimension, which is M
        working_by_age_diff_mean=mean(working_by_age_diff,dims=2)
        working_by_experience_diff_mean=mean(working_by_experience_diff,dims=2)
        wages_by_educ_diff_mean=mean(wages_by_educ_diff,dims=2)
        wages_by_age_diff_mean=mean(wages_by_age_diff,dims=2)
        transition_diff_mean=mean(transition_diff,dims=3)

        return -working_by_educ_diff_mean, -working_by_age_diff_mean, -working_by_experience_diff_mean, -wages_by_educ_diff_mean, -wages_by_age_diff_mean, -transition_diff_mean
end





