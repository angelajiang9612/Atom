
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

df_first=filter(row -> row.age==45, df)

maximum(df_first[:,4])

N=1506
A_final=64 ##A for terminating period
A_init=45
T=20
M=20 ##S is number of simulations
#set initial parameters
##sigma eta is randomly set (=rho), b is randomly set, the rest is using results from Wholpin paper
δ=0.95

α_1=-1375
α_2=-0.047
α_3=-11.5
α_4=-95.6
b=0.1

β_1=-0.280
β_2=0.024
β_3=-0.0002
β_4=0.050


σ_η=0.038
σ_ξ=0.194

κ=100


##Grid for state variable work experience

K=45 ##the maximum first period experience level, which is 24, +20, K is length so add 1
k_grid=Int.(collect(range(0, stop = 44, length = K)))


##generate initial set of shocks

Random.seed!(1234);

xi_matrix = reshape(rand(Normal(0, σ_ξ), N*T*M), N,T,M)

xi_matrix

####

########compute cutoffs#########

## Calculate the required moments from the data

#1 Average periods working by education

df_data=filter(row -> row.age<=54, df)
df_data_edu_under_12=filter(row -> row.age<=54 && row.educ<=11, df)
df_data_edu_12=filter(row -> row.age<=54 && row.educ==12, df)
df_data_edu_13_15=filter(row -> row.age<=54 && 13<=row.educ<=15, df)
df_data_edu_16_plus=filter(row -> row.age<=54 && row.educ>=16, df)

Working_Average=mean(df_data[:,3]) #overall working around 59% of the time
Working_Average_under_12=mean(df_data_edu_under_12.lfp) #overall working around 59% of the time
Working_Average_12=mean(df_data_edu_12.lfp) #overall working around 59% of the time
Working_Average_13_15=mean(df_data_edu_13_15.lfp)
Working_Average_16_plus=mean(df_data_edu_16_plus.lfp)


#2 Average periods working by age

function working_by_age()
    participation_by_age=zeros(10)
    for t=45:54
        df_age=filter(row -> row.age==t, df_data)
        participation_by_age[t-44]=mean(df_age.lfp)
    end
    return participation_by_age
end

participation_by_age=working_by_age()

#in this dataset, participation by age is the same for all the 10 years of data we have

#3 Fraction of working by work experience

df_data_exp_10_or_less=filter(row -> row.x<=10, df_data)
df_data_exp_11_to_20=filter(row -> 11<=row.x<=20, df_data)
df_data_exp_21_plus=filter(row -> row.x >= 21, df_data)

Working_exp_10_or_less=mean(df_data_exp_10_or_less.lfp)
Working_exp_11_to_20=mean(df_data_exp_11_to_20.lfp)
Working_exp_21_plus=mean(df_data_exp_21_plus.lfp)

#4 Mean wage by categories


df_data
df_data_working=filter(row -> row.lfp==1, df_data)
df_data_working_educ_under_12=filter(row -> row.educ<=12, df_data_working)
df_data_working_educ_12=filter(row -> 11<=row.educ==12, df_data_working)
df_data_working_educ_13_15=filter(row -> 13 <= row.educ <= 15,  df_data_working)
df_data_working_educ_16_plus=filter(row -> row.educ>=16,  df_data_working)

Wages_average=mean(df_data_working.wage)
Wage_Average_educ_under_12=mean(df_data_working_educ_under_12.wage) #overall working around 59% of the time
Wage_Average_educ_12=mean(df_data_working_educ_12.wage) #overall working around 59% of the time
Wage_Average_educ_13_15=mean(df_data_working_educ_13_15.wage)
Wage_Average_educ_16_plus=mean(df_data_working_educ_16_plus.wage)


function wages_by_age()
    wages_by_age=zeros(10)
    for t=45:54
        df_age=filter(row -> row.age==t, df_data_working)
        wages_by_age[t-44]=mean(df_age.wage)
    end
    return wages_by_age
end

wages_age=wages_by_age()




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

Transition_Matrix



function Transition_Estimate(df) #this function calculates transition matrix
    for l=1:2

    end
    return counter_0, counter_1
end

a,b=Transition_Estimate(df)


a









##

#Data Moments Results

#1.

Working_Average
Wage_Average_educ_under_12
Working_Average_13_15
Working_Average_16_plus

#2 Average periods working by age

participation_by_age

#3 Fraction of working by work experience

Working_exp_10_or_less
Working_exp_11_to_20
Working_exp_21_plus

#4 Mean wage overall, by schooling, by age


Wages_average

Wage_Average_educ_under_12
Wage_Average_educ_12
Wage_Average_educ_13_15
Wage_Average_educ_16_plus

wages_age


#5

Transition_Matrix


##change this

function cutoffs_analytical()
    EV=zeros(N,K,T)
    xi_star=zeros(N,K,T)  ##some of these will be zeros because they are not calculated...
    for i=1:N
        println("this is", i, "person")
        for t=T:-1:1
            xi=xi_matrix[i,t,:]
            row=Matrix(filter(row -> row.id ==i && row.age==44+t, df))
            y=row[8]
            s=row[6]
            V_T=zeros(N,M,K)
            if t==T
                for k=1:K #experience
                    k_last=k_grid[k]
                    E_V_k_T=mean(V_T[i,:,k])
                    EV[i,k,t]=E_V_k_T
                    xi_star[i,k,t]=log(-α_1 - α_2*y + b*(1+α_2) - α_3*k_last -α_4*s) - (β_1 + β_2*k_last + β_3*k_last^2 + β_4*s) - log(1+α_2)
                end
            elseif t < T
                for k=1:K-t ##need to -t, otherwise cannot always get next period.
                    k_last=k_grid[k]
                    E_V_k_t=mean(V_T[i,:,k])
                    EV[i,k,t]=E_V_k_t
                    xi_star[i,k,t]=log((-α_1 - α_2*y + b*(1+α_2) - α_3*k_last -α_4*s) + δ*(EV[i,k,t+1]-EV[i,k+1,t+1])) - (β_1 + β_2*k_last + β_3*k_last^2 + β_4*s) - log(1+α_2)
                end
            end
        end
    end
    return EV, xi_star
end



##check why the xi_stars are weird


#do this for one simulation
#xi N,T,M


function getendogenous() #this function calculates endogenous choices
    k_vector=zeros(20) #keep track of the experience throughout the 20 years
    x=zeros(20) #x is the 0 1 labor decision in the period
    for i=1:N
        row=Matrix(filter(row -> row.id ==i && row.age==45, df))
        for t=1:20
            xi=xi_matrix[i,t,1] #i's t period shock, just use the first simulation for now
            if t==1
                k_vector[t]=row[4] #initial period experience just use the one in data
                k= row[4]
                xi_cutoff=xi_star[i,k,t]
                if xi >= xi_cutoff
                    k_vector[t+1]=k+1 #update experience next period
                    x[t]=1 #this period participation is 1
                elseif xi < xi_cutoff
                    k_vector[t+1]=k
                    x[t]=0
                end
            elseif 1<t<20
                k= floor(Int,k_vector[t])
                xi_cutoff=xi_star[i,k,t]
                if xi >= xi_cutoff
                    k_vector[t+1]=k+1 #update experience next period
                    x[t]=1 #this period participation is 1
                elseif xi < xi_cutoff
                    k_vector[t+1]=k
                    x[t]=0
                end
            else t==20 #do not need to update next period.
                k= floor(Int,k_vector[t])
                xi_cutoff=xi_star[i,k,t]
                if xi >= xi_cutoff
                    x[t]=1 #this period participation is 1
                elseif xi < xi_cutoff
                    x[t]=0
                end
            end
        end
        return x, k_vector
    end
end

x,k=getendogenous()



fffff

zzkjj


using code
lalalala

ssss
