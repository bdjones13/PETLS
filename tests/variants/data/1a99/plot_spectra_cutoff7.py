import pandas as pd
df = pd.read_csv("./1a99_cutoff7/spectra_summary.txt",sep="\t")




import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [25, 15]
plt.rcParams['figure.dpi'] = 100 # for test quality, use 100, publication 600
sns.set_theme()
sns.set_style("whitegrid")
sns.color_palette("colorblind")
sns.set(rc={'text.usetex': True})

fig,axes = plt.subplots(3,2)

# graph display parameters
fontsize=22
y_lims_betti = (-0.5,100)
y_lims_lambda = (-0.5,4) # best to have them be the same if scale allows 
x_lims = (-0.1,7.1)
label_offset_x = -1.25
label_offset_y_betti = (y_lims_betti[0]+y_lims_betti[1])/2
label_offset_y_lambda = (y_lims_lambda[0]+y_lims_lambda[1])/2
yticks_betti = [10*i for i in range(11)]
yticklabels_betti = yticks_betti
yticks_lambda = [i for i in range(6)]
yticklabels_lambda = yticks_lambda


def plot_spectra(axes_row,axes_col,df,y_col,y_label,marker,y_lims,label_offset_y,yticks,yticklabels):
    sns.lineplot(ax=axes[axes_row,axes_col],data=df,x="filtration",y=y_col,color="blue",
                 drawstyle='steps-post',marker=marker,markersize=0, markerfacecolor='white',
                 markeredgecolor='blue',markeredgewidth=0,
                 )
    axes[axes_row,axes_col].text(label_offset_x,label_offset_y,y_label, ha="left",va="center",fontsize=fontsize,rotation=0)
    axes[axes_row,axes_col].set_ylabel("")
    axes[axes_row,axes_col].set_xlabel("Filtration",fontsize=fontsize)
    axes[axes_row,axes_col].set_ylim(y_lims)
    axes[axes_row,axes_col].set_xlim(x_lims)
    axes[axes_row,axes_col].set_xticks([0,1,2,3,4,5,6,7])
    axes[axes_row,axes_col].set_yticks(yticks)
    axes[axes_row,axes_col].set_xticklabels([0,1,2,3,4,5,6,7],fontsize=fontsize)
    axes[axes_row,axes_col].set_yticklabels(yticklabels,fontsize=fontsize)
    


betti_args = ['s',y_lims_betti,label_offset_y_betti,yticks_betti,yticklabels_betti]
lambda_args = ['o',y_lims_lambda,label_offset_y_lambda,yticks_lambda,yticklabels_lambda]
plot_spectra(0,0,df,"betti_0",r"$\beta_0^{a,b}$",*betti_args)    
plot_spectra(1,0,df,"betti_1",r"$\beta_1^{a,b}$",*betti_args)


y_lims_betti = (-0.5,400)
label_offset_y_betti = (y_lims_betti[0]+y_lims_betti[1])/2
label_offset_y_lambda = (y_lims_lambda[0]+y_lims_lambda[1])/2
yticks_betti = [40*i for i in range(11)]
yticklabels_betti = yticks_betti
betti_args = ['s',y_lims_betti,label_offset_y_betti,yticks_betti,yticklabels_betti]
plot_spectra(2,0,df,"betti_2",r"$\beta_2^{a,b}$",*betti_args)   

plot_spectra(0,1,df,"lambda_0",r"$\lambda_0^{a,b}$",*lambda_args)    
plot_spectra(1,1,df,"lambda_1",r"$\lambda_1^{a,b}$",*lambda_args)
plot_spectra(2,1,df,"lambda_2",r"$\lambda_2^{a,b}$",*lambda_args)        

#fill under graph
for type_id in range(0,2):
    for dim in range(0,3):        
        l = axes[dim,type_id].lines[0]
        x = l.get_xydata()[:,0]
        y = l.get_xydata()[:,1]
        axes[dim,type_id].fill_between(x,y,alpha=0.6,step="post")

plt.subplots_adjust(hspace=0.4)
plt.savefig("1a99_cutoff7_spectra.png")