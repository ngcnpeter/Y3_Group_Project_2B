from TCSPC import *

def n_case_df(df_list,col):
    '''Return df for cases in df_list
       Input:
       df_list  list of the structure [[df_1,df_2,...df_n_case],...,] (nested list of 20 (n_photon_arr) by n_case)
       col      col names for the cases'''
    df_list_case = []
    for df in df_list:
        n_val_df =pd.concat(df,keys = col,axis = 1) #concat dfs of n cases along axis 1
        df_list_case.append(n_val_df) #append each photon number case
    return pd.concat(df_list_case,keys = range(20)) 

EGFP = Phasor([0.497,0.503],[2.43,3.07]) #Initialize phasor object

#############  LOAD DATA ########################
EGFP.y_list,EGFP.y_bg_list,EGFP.phasor_list,EGFP.phasor_bg_list = np.load('df/EGFP_data.npy')
EGFP.y_list=EGFP.y_list.astype(int)
EGFP.y_bg_list=EGFP.y_bg_list.astype(int)

#############  REPEAT FITTING/PARAMETER GENERATION FOR DATA ########################
n_photon_arr = np.logspace(4,9,20)
df_list = []   #mle
df_ls_list = [] #least squares
df_p_list = [] #phasor
N=2 #N components
for i in range(len(n_photon_arr)):
    df_bg = EGFP.val_df(N,sim_data=EGFP.y_bg_list[i])
    df_no_bg = EGFP.val_df(N,sim_data=EGFP.y_list[i],bg=False)
    df_rescale= EGFP.val_df(N,sim_data=EGFP.y_list[i],rescale = True)
    df_no_bg_rescale= EGFP.val_df(N,sim_data=EGFP.y_list[i],r=30,rescale = True,bg=False)
    df_list.append([df_bg,df_no_bg,df_rescale,df_no_bg_rescale])

for i in range(len(n_photon_arr)):
    df_bg = EGFP.val_df(N,resid_func = LS_deviance_residual,sim_data=EGFP.y_bg_list[i])
    df_no_bg = EGFP.val_df(N,resid_func = LS_deviance_residual,sim_data=EGFP.y_list[i],bg=False)
    df_rescale= EGFP.val_df(N,resid_func = LS_deviance_residual,sim_data=EGFP.y_list[i],rescale = True)
    df_no_bg_rescale= EGFP.val_df(N,resid_func = LS_deviance_residual,sim_data=EGFP.y_list[i],r=30,rescale = True,bg=False)
    df_ls_list.append([df_bg,df_no_bg,df_rescale,df_no_bg_rescale])

for i in range(len(n_photon_arr)):
    df_bg_1 = EGFP.generate_df(phasor_data = EGFP.phasor_bg_list[i],idx = None)
    df_bg_2 = EGFP.generate_df(phasor_data = EGFP.phasor_bg_list[i],idx = [0,2,4]) #use harmonics aside from first
    df_no_bg = EGFP.generate_df(phasor_data = EGFP.phasor_list[i])
    df_p_list.append([df_bg_1,df_bg_2,df_no_bg])

#################  SAVE DATAFRAMES ########################
df_mle = n_case_df(df_list,col = ['bg','no_bg','bg_rescale','no_bg_rescale']) #mle
df_ls  = n_case_df(df_ls_list,col = ['bg','no_bg','bg_rescale','no_bg_rescale']) #leastsq
df_p   = n_case_df(df_p_list,col = ['bg1','bg2','no_bg'])                       #phasor

df_mle.to_csv('df/df_mle.csv')
df_ls.to_csv('df/df_ls.csv')
df_p.to_csv('df/df_p.csv')