# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:17:07 2012

@author: michelangelo puliga
@email: puligam@ethz.ch

Description: Compute the DebtRank (S.Battiston et al.) within the networkX framework for 
network analysis and visualization.

Dependencies: matplolib, networkx, colorsys

Input: network file in edgelist (weighted) format
Output: vector of the debtrank and the spiral plot of the DR


"""

import networkx as nx
import sys
import glob
import numpy as np
import math
import re

import matplotlib.pyplot as plt
from colorsys import *
from pylab import *

from matplotlib.path import Path
import matplotlib.patches as patches

from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator

class work:
    """
    main simulations and plots
    
    """
    
    
    def __init__(self,nodef):
        self.filenode = nodef  ## the file containing the variables
        
    def toDate(self, filn):
        """
        This function converts the name of the output file (i.e. 2005-2_DR.png ... 2005-2)
        in an ordinal (datetime.toordinal) number (see python datetime doccumentation)
        The output is needed by the plot_date function
        """
        filn = filn.replace("_DR.png","")
        sp = filn.split("-")
        t = datetime.datetime(int(sp[0]), int(sp[1]), 1)
        ret = t.toordinal()
        return ret
        
    def loadNetwork(self,filename, filterlow, filterup):
        """
        load the network: must be in node1, node2, weight format
        i.e.
        
        1  2  0.33
        1  3  0.35
        
        NOTE. here you can add a 4th property to the nodes an economic size
        changing with time
        """
        
        print filename
        op = open(filename,"rb")
        self.G = nx.DiGraph() ## directed graph
        for l in op.readlines():
            l = l[:-1]
            sp = l.split()
            id0 = int(sp[0].strip())
            id1 = int(sp[1].strip())
            weight = float(sp[2])
            if(weight > filterlow and weight < filterup):
                self.G.add_edge(id0,id1,impact=weight)
        op.close()
        
        ## because the adjacency list can be incomplete load 
        ## the remaining nodes from the node file
        op = open(self.filenode,"rb")
        for l in op.readlines():
            l = l[:-1]
            sp = l.split(",")
            id0 = int(sp[0].strip())
            oldnodes = self.G.nodes()
            try: # add the node only if it is not already present
                oldnodes.index(id0)
            except: 
                self.G.add_node(id0)
        op.close()
        
        

                    

    def addNodeAttributes(self,net,attributes, attribute_name):
        """
        add attributes to the node (i.e. add the "size", or others)
        """
        for k in attributes.keys():
            net.node[k][attribute_name] = attributes[k]


    def addAttributesForDebtRank(self,G,attribute_list):
        """
        add all the attributes for the debt rank (h, s, size)
        """
        for attr in attribute_list:
            self.addNodeAttributes(G,attr.values()[0], attr.keys()[0])
    

    def seedNodes(self,G,num_seeds = None, list_seeds = None, random_nodes = False):
        """
        seed the nodes to initialize the DR simulation.
        num_seeds = number of nodes to be impacted at the beginning of the simulation. It is active only for random picking of the nodes
        list_seeds = vector of the nodes to be impacted at the beginning of the simulation
        random = boolean for picking or not random nodes
        
        DEFAULT. 
        initial impact = 1
        only one node at time is impacted
        """

        ## impact the nodes one at time
        self.nodes_seed = []

        if(group_debt_rank):
            for n in G.nodes():
                G.node[n]['h'] = psi # psi is the initial impact ( 0 <psi<1)
                G.node[n]['s'] = 1
        else:
            if(random_nodes):
                nodes_tmp = np.random.randint(1,len(self.G.nodes()),num_seeds)
                for i in nodes_tmp:
                    self.nodes_seed.append(G.nodes()[i])
                for n in self.nodes_seed:
                    G.node[n]['h'] = 1. # initial impact equal to 1.
                    G.node[n]['s'] = 1  # update the status of an impacted node
            else:
                for n in G.nodes():
                    G.node[n]['h'] = 0.
                    G.node[n]['s'] = 0
                for n in list_seeds:
                    G.node[n]['h'] = 1.
                    G.node[n]['s'] = 1
                
        ##
                
        self.original_h= {}
        for n in G.nodes():
            self.original_h[n] = G.node[n]['h']
        self.active_nodes =[]
        for i in G.nodes():
            if(G.node[i]['h'] > 0 and G.node[i]['s'] < 3):
                self.active_nodes.append(i)
            

		
    def startImpactDynamics(self,G,active_nodes):
        """
        compute the forward impact
        from a node j in the active nodes
        count the impact to the other nodes
        """
        max_steps = 100
        count_steps = 0
        num_active_nodes = len(active_nodes) 
        while(num_active_nodes > 0 and count_steps < max_steps):
            for j in active_nodes:
                succ_j = G.successors(j)
                ## compute the shocks caused by j
                shocks = {}
                for suc in succ_j:
                    Wsj = G[j][suc]['impact']
                    Hp = G.node[j]['h']
                    shocks[suc] = Wsj*Hp
                ## propagate the shock                    
                for suc in succ_j:                    
                    G.node[suc]['h'] = min(1.,G.node[suc]['h']+shocks[suc])
                ## update the shock
            ### update the state of the nodes
            for suc in G.nodes():
                if(G.node[suc]['s'] == 2):
                    G.node[suc]['s'] = 3
                if(G.node[suc]['s'] == 1):
                    G.node[suc]['s'] = 2

            ### recompute the active nodes 
            tot_actives = 0
            active_nodes = []
            avg_h = 0
            for i in G.nodes():
                if(G.node[i]['s'] == 0 and G.node[i]['h'] > 0):
                    active_nodes.append(i)
                    avg_h += G.node[i]['h']
                    G.node[i]['s'] = 1
                    tot_actives += 1
            
            num_active_nodes = tot_actives
            count_steps += 1
		 

    def getDebtRank(self,G,original_h):
        """
        from the original impact matrix 
        compute the debt rank
         
        Note. Include in the following lines the temporal dependence
        with the following technique
        
        add an attribute to the network nodes (i.e a "tsize" depending on time 
        load during the network creation. )
        
        
        """
        
        self.debtrank = {}        
        nodes = G.nodes()
        

        sizes = []
        for n in nodes:
            sizes.append(G.node[n]['size'])
        
        sumsize = np.sum(sizes)
        ## normalize the sizes
        for n in nodes:
            vfin = G.node[n]['h'] * G.node[n]['size']/sumsize ## here add the economic size depending on time G.node[n]['tsize'] 
            vini = original_h[n] * G.node[n]['size']/sumsize ## here add the economic size depending on time G.node[n]['tsize']
            self.debtrank[n] = vfin - vini



    def plotSeriesOfDR(self,drseries,dates,numnodes):
        """
        Plot the DebtRank over time for all the nodes with color ordered by DR value
        
        input: dr series, dates, and numof nodes
        """

        def torgb(h,s,v):
            # convert hsv triplet into rgb codes
            (r,g,b) = hsv_to_rgb(h,s,v)
            rgb = (r,g,b)
            return rgb        
    
    
        def colorRainbow(numcolors):
            # create the rgb triplets 
            vals = np.linspace(0.1,0.6666,numcolors)
            triplets = []
            for v in vals:
                triplets.append(torgb(v,1.,1.))        
            return triplets       
        
        
        
        
        yaxis = "Debt Rank"
        title = "EMID DR"
        tseries = {}


        colors = colorRainbow(numnodes)
        
            
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_ylabel(yaxis, fontsize=14)

#        ax.set_autoscale_on(False)        
        ax.set_ylim(0., 1.)
#        ax.set_ylim(1, 150)
#        ax.set_xscale('log')
#        ax.set_yscale('log')

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)


        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)



        i = 0
        for lab in drseries.keys():
            
            try:
                ax.plot_date(dates,drseries[i],  color=colors[i],  linestyle="-" , linewidth=2,  alpha=0.7, marker = None)
            except:
                print len(dates), len(drseries[i])
            i += 1
        ax.xaxis.set_major_locator( mdt.MonthLocator(interval = 20))
        ax.xaxis.set_major_formatter(mdt.DateFormatter('%b\n%Y'))
        #handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels, legend_position,ncol=6)
        #ax.legend(handles, labels, 'right', bbox_to_anchor=(1.5, 1.05),ncol=2, fancybox=True, shadow=True)
        
        plt.title(title,  fontsize = 14)

        
        plt.grid()         
                               
            
    def plotDR(self,G,title, fileout = None):
        """
        Plot the spiral of the debt rank, input the NetworkX graph (h,s,size,label) additional node parameters)
        """
        def torgb(h,s,v):
            """
            transform from linear colorspace hsv to rgb palette
            """
            (r,g,b) = hsv_to_rgb(h,s,v)
            rgb = (r,g,b)
            return rgb        
    
    
        def colorRainbow(numcolors):
            """
            create the rgb triplets for the plot
            """
            vals = np.linspace(0.1,0.6666,numcolors) ## changing the values here change the colorscale (use it careful, it is hsv space)
            triplets = []
            for v in vals:
                triplets.append(torgb(v,1.,1.))        
            return triplets
        
        def spiral(laps):
            """
            create the spiral points
            """
            n = 1000
            spiralc = []
            for i in np.arange(n):
                r = float(i)/float(n)
                phi = float(laps) * float(r) * float(math.pi) 
                x = (1 - r) * np.cos(phi)
                y = (1 - r) * np.sin(phi)
            
                spiralc.append([x,y])
            return spiralc
        
        
        def getCoord(G,nodeid):
            """
            coordinates in spiral representation
            """
            r = G.node[nodeid]['h']
            phi = r * laps * math.pi
            xr = (1- r) * math.cos(phi)
            yr = (1- r) * math.sin(phi)
            f = (xr,yr)
            return f

            
        dranks = []
        sizes = []
        labels = []
        x = []
        y = []
        colors = []
        
        ## parameters of the plot
        laps = 20.
        numcolors = 100
        multiplier = 2.5
        ##

        for n in G.nodes():
            f = getCoord(G,n)
            x.append(f[0])
            y.append(f[1])
            dranks.append(G.node[n]['h'])
            sizes.append(float(G.node[n]['size']) * multiplier)
            labels.append(G.node[n]['label'])

        maxdr = max(dranks)
        colorsrb = colorRainbow(numcolors)
        
        lx = len(G.nodes())
        colornode = {}
        i = 0
        for dr in dranks:
            idx = int( numcolors * dr)  
            colors.append(colorsrb[numcolors - idx - 1])
            colornode[i] = colorsrb[numcolors - idx - 1]
            i += 1
        spiralc = spiral(laps)
        spiralx = []
        spiraly = []
        for c in spiralc:
            spiralx.append(c[0])
            spiraly.append(c[1])
                
        fig = plt.figure()
        ax = fig.add_subplot(111)

      
        #ax.set_ylabel(yaxis, fontsize=16)

        #ax.set_autoscale_on(False)        
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        #ax.set_xscale('log')
        #ax.set_yscale('log')

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(16)


        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(16)


        for e in G.edges():
            node1 = e[0]
            c1 = getCoord(G,node1)
            node1x = c1[0]
            node1y = c1[1]
            node2 = e[1]
            c2 = getCoord(G,node2)
            node2x = c2[0]
            node2y = c2[1]
            weight = G[node1][node2]['impact']
            wd = weight * 10.
            
            # use it to generate curved links (comment Line2D in this case) pay attention to disconnected nodes ! (out lines)
            #c = patches.FancyArrowPatch((node1x,node1y), (node2x, node2y), lod = True, arrowstyle="->", mutation_scale=wd ,connectionstyle = "arc3,rad=0.4",linewidth=1,alpha=0.4,color = colornode[node2] ,zorder=0)
            #ax.add_patch(c)
            ##
            line = Line2D((node1x,node2x), (node1y, node2y),alpha=0.1, lw=0.5,color = colornode[node2] ,zorder=0, solid_capstyle='butt')
            ax.add_line(line)
            

        ax.plot(spiralx,spiraly, c='black',linestyle=":" , linewidth=1)
        ax.scatter(x,y,sizes, c=colors,zorder=3)
        for n in np.arange(len(G.nodes())):
            ax.annotate(labels[n],(x[n],y[n]), fontsize = 8)
        
        plt.title(title,  fontsize = 16) ## add here the title
        
        #plt.grid() 
        plt.savefig(fileout,dpi=200, figsize=(12, 8), bbox_inches='tight',pad_inches=0.3) ## save the picture to png
        plt.close(fig) ## remember to close the figure
        


def getSize(nodeid, timevariable = None):
    """
    here add the code for changing the size 
    eventually during the simulation
    """
    size = 20
    return size
    
def getLabelSize(filenode):
    """
    here add the code for getting the label 
    in this case the name of the country
    format: index,countryname
        
    MA major
    GR large
    ME medium
    PI small
    MI minor
    
    Where "Banks are classified into 5 groups according to their weighted
    asset portfolio: major banks (higher than 60 billion euro), large
    banks (from 26 to 60 billion euro), medium banks (from 9 to 26 billion
    euro), small banks (from 1.3 to 9 billion euro) and minor banks (less
    than 1.3 billion euro)."
    
    """
    
    sizes = {'MA' : 80., 'GR' : 43., 'ME' : 17. ,'PI': 5. , 'MI' : 0.5, 'ND' : 17.} ## average.
    
    op = open(filenode,'r')
    labelsize  = []
    for o in op.readlines():
        o = o[:-1]
        sp = o.split(",")
        labelsize.append([sp[1],sizes[sp[2].strip()]])
    op.close()
    return labelsize
    

def inizializeNodes(wk):
    """
    function to initialize the networkx graph with the attributes 
    needed by the DR simulation:
    
    h, s (impact and status)
    size (size of the node)
    label (label attached to the node)
    """
    
    attribute_list = []
    
    attribute_name = "h"
    attributes = {}
    for i in wk.G.nodes():
        attributes[i] = 0.
    
    attribute_list.append({attribute_name : attributes})
        
    attribute_name = "s"
    attributes = {}
    for i in wk.G.nodes():
        attributes[i] = 0
    
    attribute_list.append({attribute_name : attributes})
    
    labelsize = getLabelSize(wk.filenode)
    attribute_name = "label"
    attributes = {}
    

    for i in wk.G.nodes():
        attributes[i] = labelsize[i][0]
    
    attribute_list.append({attribute_name : attributes})
    

    attribute_name = "size"
    attributes = {}
    
    for i in wk.G.nodes():
        attributes[i] = labelsize[i][1]
    
    attribute_list.append({attribute_name : attributes})
    
        
    return attribute_list
    


def computeDr(wk,netfile, filterlow = -100., filterup = 100.):  
    """
    compute the Debt Rank (impacting the nodes one at time with initial impact = 1)
    filter the network according to filterlow and filterup limits (default -100, 100, quite large limits)
    """

    net = netfile    
    wk.loadNetwork(net,filterlow,filterup) ## load the network with filtering
    attributes = inizializeNodes(wk)
    
    wk.addAttributesForDebtRank(wk.G,attributes) # add the attributes for the DR simulation to the graph
    
    drank = {}
    if (group_debt_rank):
        wk.seedNodes(wk.G,1,[node],False) # seed nodes once by psi

    ## run through all nodes impacting one node at time with h = 1
    for node in wk.G.nodes():
        if(not group_debt_rank):
            wk.seedNodes(wk.G,1,[node],False) # seed nodes (one at time)
        wk.startImpactDynamics(wk.G,wk.active_nodes) # impact dynamics
        wk.getDebtRank(wk.G,wk.original_h) # debt rank (normalized impact)
        drank[node] = np.sum(wk.debtrank.values()) # debt rank as the total size of the effects due to a single node
    return drank





########################################
########### MAIN #######################
########################################

# start here the simulation

#initialize the class for work
wk = work("bank-info.csv") ## pass it the node file


# get the debt rank

gfiles = glob.glob("graphs/*eMid") ## use the glob module and function to fetch the network files

ddates = []
nnodes = 313

## in case of group debt rank use this (set it to false otherwise)
group_debt_rank = True
psi = 0.1

drseries = {}

for gf in sorted(gfiles):
    fileout = gf.replace("-weighted-edgelist-aggregated-eMid","_DR.png")
    date = wk.toDate(fileout.replace("graphs/",""))
    ddates.append(date)

for n in np.arange(nnodes):
    drseries[n] = []



for gf in sorted(gfiles):
    drank = computeDr(wk,gf) ## change here the file for network analysis
    
    fileout = gf.replace("-weighted-edgelist-aggregated-eMid","_DR.png") ## create the output file 
    title = gf.replace("-weighted-edgelist-aggregated-eMid"," Debtrank for emid")
    title = title.replace("graphs/","")
    

    
    for node in drank.keys():
        #print "Node", node, "Debt Rank", drank[node] ### add here the output 
        wk.G.node[node]['h'] = drank[node] ## note at the *end* of the simulation for simplicity we replaced the impact "h" with the debt rank
        drseries[node].append(drank[node])
        
    #### call this function to visualize the debt rank
    #wk.plotDR(wk.G,title,fileout)

### plot the series of the debt rank
wk.plotSeriesOfDR(drseries,ddates,nnodes)
    
# do not forget to show the plot
plt.show()

