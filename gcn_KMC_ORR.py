import math
import numpy as np
import collections as co
import time
import matplotlib.pyplot as plt
import os


#to_read= raw_input("Insert the name of the file containing the catalyst you wish to use for this simulation: ")
#type(to_read)
to_read= "Geometrical.xyz"
start = time.time()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
    
length=file_len(to_read)

counterline=0
countercluster=0

os.system("mkdir Output")
with open("gcn_bridge_genome.dat", "w")as o:
	with open(to_read, "r") as f:
		while counterline <length:
			countercluster+=1
			natoms=int(f.readline())
			secondline=f.readline()
			x=[]
			y=[]
			z=[]
			numN=[]
			numN=[0]*natoms
			lines=[]
			cutoffNN=3.88*math.sqrt(2)/2*1.2
			cn_max=18.0
			num_pair=0
			for i in range(natoms):
				lines.append(f.readline())
			counterline+=2
			for i in lines:
				counterline+=1
				x.append(float(i.split()[1]))
				y.append(float(i.split()[2]))
				z.append(float(i.split()[3]))
				
				
			def dist(i,j):
				return math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
			
			
			
			
			for i in range(natoms-1):
				for j in range(i+1, natoms):			
					if dist(i,j)<cutoffNN:
						numN[i]+=1
						numN[j]+=1
			def gcn():
				gcn=[]
				
				for i in range(natoms-1):
					for j in range(i+1,natoms):
						if dist(i,j)<cutoffNN:
							neighb=[]
							for k in range(natoms):
								if i!=k and j!=k:
									if dist(i,k)<cutoffNN or dist(j,k)<cutoffNN:
										neighb.append(numN[k])
							neighbours=np.asarray(neighb)			
							if numN[i] <10 and numN[j]<10:
								gcn.append(np.sum(neighbours)/cn_max)
							else:
								gcn.append(0)
			
				return gcn
			
			
			
			GCN=gcn()	
			
			
		
			
			counter=co.Counter(GCN)
			
			
			
			def cg1():
				nl=0
				nm=0
				nh=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl+=counter[key]
						if key>5.5 and key<7.5:
							nm+=counter[key]
						if key>=7.5:
							nh+=counter[key]
				
				cg_list=[nl, nm, nh]
				return cg_list
							
				
				
			def cg2():
				nl=0
				nm=0
				nh=0
				nsh=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl+=counter[key]
						elif key>=8.3:
							nsh+=counter[key]
						elif (key>=7.5 and key<8.3) or (key>=6.59 and key<=6.61):
							nh+=counter[key]
						else:
							nm+=counter[key]
				
				cg_list=[nl, nm, nh, nsh]
				return cg_list
				
      def cg3():
				nl=0
        nlm=0
				nmm=0
        nmh=0
				nhh=0
				for key in counter:
					if key!=0:
						if key<=5.2:
							nl+=counter[key]
						if key>5.2 and key<6:
							nlm+=counter[key]
            if (key>6 and key<6.6) or (key>6.7 and key<7):
              nmm+=counter[key]
            if (key>6.6 and key<6.63) or (key>7.1 and key<7.3):
              nmh+=counter[key]
						if (key>6.63 and key<6.7) or (key>7.3):
							nhh+=counter[key]
			        nm = nlm + nmm
              nh = nmh + nhh
				cg_list=[nl, nm, nh]
				return cg_list
				
				
			def cg1_and_2():
				nl1=0
				nm1=0
				nh1=0
				nsh1=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl1+=counter[key]
						if key>5.5 and key<7.5:
							nm1+=counter[key]
						if key>=7.5:
							nh1+=counter[key]
							
				nl2=0
				nm2=0
				nh2=0
				nsh2=0
				for key in counter:
					if key!=0:
						if key<=5.5:
							nl2+=counter[key]
						elif key>=8.3:
							nsh2+=counter[key]
						elif (key>=7.5 and key<8.3) or (key>=6.59 and key<=6.61):
							nh2+=counter[key]
						else:
							nm2+=counter[key]		
				
				print "Cluster #",countercluster
				print "\nLow gcn sites: ", nl1, "  ",nl2
				print "\nMedium gcn sites: ", nm1, "  ",nm2
				print "\nHigh gcn sites: ", nh1, "  ",nh2
				print "\nSuper-high gcn sites: ", nsh2
				print "\nTotal surface sites: ", nl1+nm1+nh1+nsh1
				print "\n**********************************************************"
				
				if countercluster==1:
					o.write("Coarse graining mode=1 and 2")
					
				
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl1)+", "+str(nl2))
				o.write("\n\nMedium gcn sites: "+ str( nm1)+", "+str(nm2))
				o.write("\n\nHigh gcn sites: "+ str( nh1)+", "+str(nh2))
				o.write("\n\nSuper-high gcn sites: "+ str( nsh1))
				o.write("\n\nTotal surface sites: "+ str( nl1+nm1+nh1))
				o.write("\n\n********************************************************")
							
				#y1=[nl1,nm1,nh1,nsh1]
				y1=[nl1*100.0/(nl1+nm1+nh1+nsh1),nm1*100.0/(nl1+nm1+nh1+nsh1),nh1*100.0/(nl1+nm1+nh1+nsh1),nsh1*100.0/(nl1+nm1+nh1+nsh1)]
				#y2=[nl2,nm2,nh2,nsh2]
				y2=[nl2*100.0/(nl2+nm2+nh2+nsh2),nm2*100.0/(nl2+nm2+nh2+nsh2),nh2*100.0/(nl2+nm2+nh2+nsh2),nsh2*100.0/(nl2+nm2+nh2+nsh2)]
				x=["Low gcn","Medium gcn","High gcn","Super-high gcn"]
				y_pos=np.arange(4)
				plt.bar(y_pos, y1, color="b", width=0.25 , align="edge")
				plt.bar(y_pos+0.25, y2, color="r", width=0.25, align="edge")
				plt.ylabel("Percentage of sites")
				#plt.ylabel("Number of sites")
				plt.xticks(y_pos, x)
				plt.title("GCN coarse-graining")
				#plt.text("Low gcn sites: ", nl1," ", nl2,"\nMedium gcn sites: ", nm1," ", nm2,"\nHigh gcn sites: ", nh1," ", nh2,"\nSuper-high gcn sites: ", nsh1)
				plt.savefig("snap"+str(countercluster))
				plt.close()			
			
			
			#Which coarse-graining?
			cg=3
			
			if cg==1:
				coarse_grain=cg1()
				nl=coarse_grain[0]
				nm=coarse_grain[1]
				nh=coarse_grain[2]
				print "Cluster #",countercluster
				print "\nLow gcn sites: ", nl
				print "\nMedium gcn sites: ", nm
				print "\nHigh gcn sites: ", nh
				print "\nTotal surface sites: ", nl+nm+nh
				print "\n**********************************************************"
				
				if countercluster==1:
					o.write("Coarse graining mode=1")
				
				with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nm)+"\n")
					p.write(str(nh)+"\n")
					p.write(str(countercluster))
				
	
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl))
				o.write("\n\nMedium gcn sites: "+ str( nm))
				o.write("\n\nHigh gcn sites: "+ str( nh))
				o.write("\n\nTotal surface sites: "+ str( nl+nm+nh))
				o.write("\n\n********************************************************")
				
				
				x=["Low gcn","Medium gcn","High gcn"]
				y=[nl*100.0/(nl+nm+nh),nm*100.0/(nl+nm+nh),nh*100.0/(nl+nm+nh)]
				#y=[nl,nm,nh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="b",width=0.50,align="center")
				plt.title("GCN coarse-graining")
				plt.ylabel("Percentage of sites")
				#plt.ylabel("Number of sites")
	
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()
			
			
			elif cg==2:
				coarse_grain=cg2()
				nl=coarse_grain[0]
				nm=coarse_grain[1]
				nh=coarse_grain[2]
				nsh=coarse_grain[3]
				print "Cluster #",countercluster			
				print "\nLow gcn sites: ", nl
				print "\nMedium gcn sites: ", nm
				print "\nHigh gcn sites: ", nh
				print "\nSuper-high gcn sites: ", nsh
				print "\nTotal surface sites: ", nl+nm+nh+nsh
				print "\n**********************************************************"
				if countercluster==1:
					o.write("Coarse graining mode=2")
				
				with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nm)+"\n")
					p.write(str(nh)+"\n")
				
				o.write("\n\nCluster #"+str(countercluster))
				o.write("\n\nGCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow gcn sites: "+ str( nl))
				o.write("\n\nMedium gcn sites: "+ str( nm))
				o.write("\n\nHigh gcn sites: "+ str( nh))
				o.write("\n\nSuper-high gcn sites: "+ str( nsh))
				o.write("\n\nTotal surface sites: "+ str( nl+nm+nh))
				o.write("\n\n********************************************************")
				
				
				x=["Low gcn","Medium gcn","High gcn","Super-high gcn"]
				#y=[nl*100.0/(nl+nm+nh),nm*100.0/(nl+nm+nh+nsh),nh*100.0/(nl+nm+nh+nsh),nsh*100.0/(nl+nm+nh+nsh)]
				y=[nl,nm,nh,nsh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="r",width=0.50, align="center")
				plt.title("GCN coarse-graining")
				#plt.ylabel("Percentage of sites")
				plt.ylabel("Number of sites")
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()
			
					
			elif cg==3:
				coarse_grain=cg3()
				nl=coarse_grain[0]
				nm=coarse_grain[1]
				#nmm=coarse_grain[2]
				#nmh=coarse_grain[3]
        nh=coarse_grain[2]
				print "Cluster #",countercluster			
				print "\nLow gcn sites: ", nl
        #print "\nMedium-low sites", nlm
				print "\nMedium-medium gcn sites: ", nm
				#print "\nMedium-high gcn sites: ", nmh
				print "\nHigh-high gcn sites: ", nh
				print "\nTotal surface sites: ", nl+nm+nh
				print "\n**********************************************************"
				if countercluster==1:
					o.write("#Coarse graining mode=3")

        with open("inputgcn.dat", "w") as p:
					p.write(str(nl)+"\n")
					p.write(str(nm)+"\n")
					p.write(str(nh)+"\n")
					p.write(str(countercluster))
				
				
				o.write("\n\n#Cluster #"+str(countercluster))
				o.write("\n\n#GCN fingerprint:\n")
				for key in counter:
				
					if key!=0:
						o.write(str("#")+str( key))
						o.write(str("	"))
						o.write(str(counter[key]))
						o.write(str("\n"))
				o.write("\n\nLow_gcn_sites = "+ str( nl))
        #o.write("\n\nLow-medium gcn sites: "+ str( nlm))
				o.write("\n\nMedium_gcn_sites = "+ str( nm))
				#o.write("\n\nMedium-high gcn sites: "+ str( nmh))
				o.write("\n\nHigh_gcn_sites = "+ str( nh))
				o.write("\n\nTotal_surface_sites = "+ str( nl+nm+nh))
				o.write("\n\n#********************************************************")
				
				
				x=["Low gcn","Medium-medium gcn","High-high gcn"]
				#y=[nl*100.0/(nl+nm+nh),nm*100.0/(nl+nm+nh+nsh),nh*100.0/(nl+nm+nh+nsh),nsh*100.0/(nl+nm+nh+nsh)]
				y=[nl,nm,nh]
				y_pos=np.arange(len(x))
				plt.bar(y_pos, y, color="r",width=0.50, align="center")
				plt.title("GCN coarse-graining")
				#plt.ylabel("Percentage of sites")
				plt.ylabel("Number of sites")
				plt.xticks(y_pos, x)
				
				plt.savefig("snap"+str(countercluster))
				plt.close()	


			elif cg==12:
				cg1_and_2()	
				
			os.system("python KMC_ORR.py")
end = time.time()
print "\n\nTime it took for the program to run"
print(end - start)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
