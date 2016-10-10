import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cPickle
#from yt.mods import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MultipleLocator
import sys


def hden_temp_hist(weight_basic,weight_col,z_key,HMkey,type=None):	
	hden_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'+z_key+'/'+HMkey+'/frbx_1kpc_500kpc_'+z_key+'_H_NumberDensity_'+weight_basic+'.cpkl'
	temp_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/'+z_key+'/'+HMkey+'/frbx_1kpc_500kpc_'+z_key+'_Temperature_'+weight_basic+'.cpkl'
	
	if type == 'column':
		column_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/'+HMkey+'/frbx_1kpc_500kpc_z02_'+weight_col+'dens.cpkl'
	elif type == 'emis':
		column_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'+z_key+'/'+HMkey+'/frbx_1kpc_500kpc_'+z_key+'_'+weight_col+'.cpkl'
	
	else:
		print "Need to give it a key so it can find the right output file!"

	hden = cPickle.load(open(hden_file,'rb'))
	temp = cPickle.load(open(temp_file,'rb'))
	coldens = cPickle.load(open(column_file,'rb'))
	
	hden,temp,coldens = np.log10(hden.flat),np.log10(temp.flat),np.log10(coldens.flat)
	
	###ORIGINAL RANGE THAT ENCOMPASSES ALL OF THE POINTS - NARROWED 10/26 FOR DIAGNOSTICS###
	#hist, xedges, yedges  = np.histogram2d(hden,temp,range=[[-6,1],[4, 6]],bins=100)#,weights=coldens)
	hist, xedges, yedges  = np.histogram2d(hden,temp,range=[[-6,-2],[4, 6]],bins=100)#,weights=coldens)	

	xbins = np.digitize(hden,xedges[1:-1])
	ybins = np.digitize(temp,yedges[1:-1])
	
	totvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
	numvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
	maxvals = np.zeros((len(xedges[1:]),len(yedges[1:])))
	
	for i in range(len(coldens)):
	    totvals[xbins[i],ybins[i]] += 10.**coldens[i]
	    numvals[xbins[i],ybins[i]] += 1
	    if 10.**coldens[i] > maxvals[xbins[i],ybins[i]]:
        	maxvals[xbins[i],ybins[i]] = 10.**coldens[i]

	avgvals = totvals/numvals
	
	idx = np.where(totvals == 0)
	avgvals[idx] = 0

	return hist,avgvals,maxvals

def add_grid(ax):
	ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
	ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
	ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
	ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
	ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
	ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75',alpha=0.2)
	ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
	ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75',alpha=0.2)
	return

try:
	plot_key = sys.argv[1]
except:
	raise Exception('hey! you need a key so it knows what plots to make! use dev or paper or evolution!')

emis = cPickle.load(open('cloudywerk.cpkl','rb'))

weight_col = 'SiIII'
weight_basic = 'SiIII_Density'
redshift_key = 'z02'
HMkey = 'g1q1'


if plot_key == 'paper':
	fileout = 'hden_temp_paper_'+redshift_key+'_'+HMkey+'ref_SiIII.pdf'
	
	#weights_col = ['SiIV','CIII','OVI']
	#weights_basic = ['SiIV_Density','CIII_Density','OVI_Density']

	weights_col = ['SiIII','SiIV','CIII']
	weights_basic = ['SiIII_Density','SiIV_Density','CIII_Density']

	#fig,ax = plt.subplots(3,2,sharex=True,sharey=True)
	#fig.set_size_inches(8,12)
	fig,ax = plt.subplots(3,1,sharex=True,sharey=True)
	fig.set_size_inches(3,9)
	ax= ax.flat	
	i,j = 0,0

	while i < len(weights_col):
		hist,avgvals,maxvals = hden_temp_hist(weights_basic[i],weights_col[i],redshift_key,HMkey,type='column')
		#im1 = ax[j].imshow(np.log10(hist.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',cmap='Blues')
		##this was j+1 but its the one that I want to keep so....
		im2 = ax[j].imshow(np.log10(avgvals.T),extent=[-6,-1,4,6],interpolation='nearest',origin='lower',vmin=9.,vmax=16.,cmap='jet')
	
		add_grid(ax[j])
		#add_grid(ax[j+1])
		if weights_col[i] != 'OVI':
			ax[j].plot(emis['hden'],emis['T'],'s',color='#6D6968',alpha=0.75,ms=3.5) ##old is HotPink
		

                ax[j].set_xlim(-6,-1)
                ax[j].set_ylim(3.5,6.)
		x0,x1 = ax[j].get_xlim()
        	y0,y1 = ax[j].get_ylim()
		ax[j].set_title(weights_col[i])

		#ax[j+1].set_aspect((x1-x0)/(y1-y0))
                #ax[j+1].set_xlim(-6,-1)
                #ax[j+1].set_ylim(3.5,6)
                #ax[j+1].set_title(weights_col[i])

	
		ax[j].set_aspect((x1-x0)/(y1-y0))
#		ax[j].set_autoscale_on(False)	
                #ax[j].set_xlim(-6,-1)
                #ax[j].set_ylim(3,6.)

		i = i + 1
		#j = j + 2
		j = j + 1

	#fig.subplots_adjust(bottom=0.125)#8)
        #cbar_ax = fig.add_axes([0.12, 0.08, 0.8, 0.01])
        #fig.colorbar(im2, cax=cbar_ax,orientation='horizontal',ticks=[10,12,14,16])#,shrink=0.6,ticks=MultipleLocator(0.5), format="%.2f")

	#for k in range(3):
	#	ax[k].set_xlim(-6.,-1.)
	#	ax[k].set_ylim(3,6.)

	#plt.ylim(3,6)
        plt.savefig(fileout)
        plt.close()

if plot_key == 'evolution':
	redshift_keys = ['z0','z02','z05','z1']
	HMkey = 'g1q1'
	
	fileout1 = 'hden_temp_evolution_emissivity_'+HMkey+'_number.pdf'
	fileout2 = 'hden_temp_evolution_emissivity_'+HMkey+'_avg.png'
	fileout3 = 'hden_temp_evolution_emissivity_'+HMkey+'_max.png'	

	weights_col = ['SiIV','CIII_977','CIV','OVI']
	weights_basic = ['Emission_SiIV','Emission_CIII_977','Emission_CIV','Emission_OVI']
	vmin,vmax = -5,5

	fig1,ax1 = plt.subplots(len(weights_col),len(redshift_keys),sharex=True,sharey=True)
	fig1.set_size_inches(16,16)
	ax1 = ax1.flat
	fig2,ax2 = plt.subplots(len(weights_col),len(redshift_keys),sharex=True,sharey=True)
        fig2.set_size_inches(16,16)
        ax2 = ax2.flat
	fig3,ax3 = plt.subplots(len(weights_col),len(redshift_keys),sharex=True,sharey=True)
        fig3.set_size_inches(16,16)
        ax3 = ax3.flat

	i,j = 0,0

	for weight in weights_col:
		i = 0
		while i < len(redshift_keys):
			print j,weight,redshift_keys[i]
			hist,avgvals,maxvals = hden_temp_hist('Emission_'+weight,weight,redshift_keys[i],HMkey,type='emis')
			im1 = ax1[j].imshow(np.log10(hist.T),extent=[-6,-2,4,6],interpolation='nearest',origin='lower',cmap='Blues')
			im2 = ax2[j].imshow(np.log10(avgvals.T),extent=[-6,-2,4,6],interpolation='nearest',origin='lower',vmin=vmin,vmax=vmax,cmap='jet')
			im3 = ax3[j].imshow(np.log10(maxvals.T),extent=[-6,-2,4,6],interpolation='nearest',origin='lower',vmin=vmin,vmax=vmax,cmap='jet')
		
			x0,x1 = ax1[j].get_xlim()
			y0,y1 = ax1[j].get_ylim()
			ax1[j].set_aspect((x1-x0)/(y1-y0))
			ax2[j].set_aspect((x1-x0)/(y1-y0))
			ax3[j].set_aspect((x1-x0)/(y1-y0))
			ax1[j].set_xticks([-6,-5,-4,-3,-2])
			ax2[j].set_xticks([-6,-5,-4,-3,-2])
			ax3[j].set_xticks([-6,-5,-4,-3,-2])

			add_grid(ax1[j])
                	add_grid(ax2[j])
			add_grid(ax3[j])	

			i = i + 1
			j = j + 1

	
	fig1.subplots_adjust(right=0.8)

        fig2.subplots_adjust(right=0.8)
        cbar_ax = fig2.add_axes([0.85,0.15,0.05,0.7])
        fig2.colorbar(im2,cax=cbar_ax)

        fig3.subplots_adjust(right=0.8)
        cbar_ax = fig3.add_axes([0.85,0.15,0.05,0.7])
        fig3.colorbar(im3,cax=cbar_ax)

	for j in range(len(redshift_keys)):
		ax1[j].set_title(redshift_keys[j])
		ax2[j].set_title(redshift_keys[j])
		ax3[j].set_title(redshift_keys[j])

	for k in range(len(weights_col)):
		ax1[k*len(redshift_keys)+len(redshift_keys)-1].text(-1,5.75,weights_col[k])
		ax2[k*len(redshift_keys)+len(redshift_keys)-1].text(-1,5.75,weights_col[k])
		ax3[k*len(redshift_keys)+len(redshift_keys)-1].text(-1,5.75,weights_col[k])
	

	fig1.savefig(fileout1)
        fig2.savefig(fileout2)
        fig3.savefig(fileout3)	

	plt.close()

		
		

if plot_key == 'dev':
	fileout = 'histdens_z02_g1q01_'+weight_basic+'_nohistweight.png'
	
	fig,ax = plt.subplots(1,3)
	fig.set_size_inches(12,4)
	ax = ax.flat

	hist,avgvals,maxvals = hden_temp_hist(weight_basic,weight_col,redshift_key,HMkey)	
	
	im1 = ax[0].imshow(np.log10(hist.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',cmap='Blues')
	im2 = ax[1].imshow(np.log10(avgvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower')#,vmin=9.,vmax=16.)
	im3 = ax[2].imshow(np.log10(maxvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower')#,vmin=9.,vmax=16.)
	
	x0,x1 = ax[2].get_xlim()
	y0,y1 = ax[2].get_ylim()
	titles = ['Number','Avg Column Density','Max Column Density']

	for j in range(len(ax)):
	   ax[j].set_aspect((x1-x0)/(y1-y0))
	   ax[j].set_xlim(-6,1)
	   ax[j].set_ylim(3.5,6)
	   ax[j].set_title(titles[j])
  
	ax[0].plot(emis['hden'],emis['T'],'s',color='HotPink',alpha=0.75,ms=3.5)

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(im3, cax=cbar_ax)#,shrink=0.6,ticks=MultipleLocator(0.5), format="%.2f")

	plt.savefig(fileout)
	plt.close()

