#!/usr/bin/env python
# coding: utf8
# copyright EDF
# author : Thierry Couvant
# Objet : modèle générique d'oxydation intergranulaire pour l'acier inoxydable
#         irradié. l'oxydation est couplée au RIS.


from math import exp, log, cos, pi
from numpy import zeros, append, isnan, mean
from random import random, sample
from pylab import figure, plot, xlabel, ylabel, title, text, legend, show, clf, grid, semilogx, axis

class LocalInitiationLaw():
    """
    Classe permettant la prévision du temps d'amorçage de la CSC des alliages de Ni par le modèle local
    """
  
    def RIS(self,a,b,c,gb,fluence,Cr_0):
        # gb[p][0]: teneur en Cr (%) à la position p du jdg, pour 0 dpa
        # gb[p][1]: teneur en Cr (%) à la position p du jdg, pour fluence
        dCr = (a + b * exp(-c * fluence))
        new_gb = []
        for p in range(len(gb)):
            new_gb.append([gb[p][0],gb[p][0]*(1+dCr)])
            
        return new_gb
    
    def PlotPureOxidationKinetics(self,tab, language):
        gbc = self.GBC
        deple = self.zCr
        
        figure(1)

        # 1: haut à droite
        # 2: haut à gauche
        # 3: bas à gauche
        # 4: bas à droite
        # 'best': optimisé par Python
        legend_location = 2
        
        clf
        if deple>0:
            semilogx(tab[1][0], tab[1][1],"r-",linewidth=3)
        semilogx(tab[0][0], tab[0][1],"b-",linewidth=3)
        if gbc>0:
            semilogx(tab[2][0], tab[2][1],"g-",linewidth=3)
 
        txt_r = str(self.cr_depleted)+"%Cr"
        txt_b = str(self.cr_nominal)+"%Cr"
        txt_g = str(self.cr_carbide)+"%Cr"
        
        if language =='english':
            xlabel(r'Time (h)')
            ylabel(r'Depth (nm)')
            if gbc>0 and deple>0:
                legend((txt_r,txt_b,txt_g),loc=legend_location)
            elif gbc>0 and deple==0:
                legend((txt_b,txt_g),loc=legend_location)
            else:    
                legend((txt_b,txt_b),loc=legend_location)
            title("Predicted intergranular oxidation kinetics of pure phases")

        if language == 'french':
            xlabel(r'Temps (h)')
            ylabel(r'Profondeur (nm)')
            if gbc>0:
                legend((txt_r,txt_b,txt_g),loc=legend_location)
            else:
                legend((txt_r,txt_b),loc=legend_location)              
            title("Cinetique prevue d'oxydation intergranulaire des phases pures")

        grid(True)
        show()

        return

    def PlotOxidationKinetics(self,pure,tab, language):

        gbc = self.GBC
        deple = self.zCr
        
        # 1: haut à droite
        # 2: haut à gauche
        # 3: bas à gauche
        # 4: bas à droite
        # 'best': optimisé par Python
        legend_location = 4

        txt_g = str(self.cr_carbide)+"%Cr"
        txt_b = str(self.cr_nominal)+"%Cr"
        txt_r = str(self.cr_depleted)+"%Cr"
        
        if deple>0:
            plot(pure[1][0], pure[1][1],"r-",linewidth=2)
        plot(pure[0][0], pure[0][1],"b-",linewidth=2)
        if gbc>0:
            plot(pure[2][0], pure[2][1],"g-",linewidth=2)

        # Cinétiques des jdg
        max_pox=0
        for gb in range(len(tab)):
            plot(tab[gb][0], tab[gb][1],"k-")
            max_pox=max(max_pox, tab[gb][1][-1])
        
        if language =='english':
            xlabel(r'Time (h)')
            ylabel(r'Predicted intergranular oxidation depth (nm)')
            if gbc>0 and deple>0:
                legend((txt_r,txt_b,txt_g,'Simulated GBs'),loc=legend_location)
            elif gbc>0 and deple==0:
                legend((txt_b,txt_g,'Simulated GBs'),loc=legend_location)
            else:
                legend((txt_b,'Simulated GBs'),loc=legend_location)
            title("Predicted intergranular oxidation kinetics")

        if language == 'french':
            xlabel(r'Temps (h)')
            ylabel(r'Profondeur (nm)')
            if gbc>0:
                legend((txt_r,txt_b,txt_g,'JdG simulés'),loc=legend_location)
            else:
                legend((txt_r,txt_b,'JdG simulés'),loc=legend_location)
            title("Cinetique prevue d'oxydation intergranulaire")

        axis([0, self.stop_time,0,max_pox])
        grid(True)
        show()

        return

    def PlotRIS(self, language):

        # 1: haut à droite
        # 2: haut à gauche
        # 3: bas à gauche
        # 4: bas à droite
        # 'best': optimisé par Python
        legend_location = 1

        txt_b = str(self.cr_nominal)+"%Cr"

        a=self.r1
        b=self.r2
        c=self.r3
        Cr_0=self.cr_nominal
        maxfluence=self.fluence
        if maxfluence==0:
            maxfluence=100

        dpa=[]
        gb0=[]

        for myfluence in range(int(maxfluence)+1):
            dpa.append(myfluence)
            dCr = a + b * exp(-c * myfluence)
            gb0.append(Cr_0 *(1+dCr))

        plot(dpa, gb0,"b-",linewidth=2)
        
        if language =='english':
            xlabel(r'Fluence (dpa)')
            ylabel(r'Intergranular Cr content (%)')
            title("Predicted RIS-induced intergranular Cr depletion")

        if language == 'french':
            xlabel(r'Fluence (dpa)')
            ylabel(r'Teneur en Cr aux joints de grains (%)')
            title("Appauvrissement en Cr induit par l'irradiation")

        grid(True)
        show()

        return
 
    def PrintGbs(self,pos,cr_nominal,cr_depleted,cr_carbide):
      
        di = []
        dg = []
        ci = []
        cg = []
        ni = []
        ng = []

        for my_gb in range(len(pos)):
            for i in range(len(pos[my_gb])):
                if pos[my_gb][i][0] == cr_carbide:
                    ci.append(i)
                    cg.append(my_gb)
                elif pos[my_gb][i][0] == cr_depleted:
                    di.append(i)
                    dg.append(my_gb)
                else:
                    ni.append(i)
                    ng.append(my_gb)

        figure()
        clf()
        plot(ng,ni,"b.")
        plot(cg,ci,"g.")
        plot(dg,di,"r.")
            
        # 1: haut à droite
        # 2: haut à gauche
        # 3: bas à gauche
        # 4: bas à droite
        # 'best': optimisé par Python
        legend_location = 4

        txt_r = str(self.cr_depleted)+"%Cr"
        txt_b = str(self.cr_nominal)+"%Cr"
        txt_g = str(self.cr_carbide)+"%Cr"
        xlabel(r'Grain boundary')
        ylabel(r'Distance from surface (nm)')
        title(r'Grain boundary set')
        legend((txt_g,txt_b,txt_r),loc=legend_location)
        show()

        
    def pureGB(self,cr_nominal,cr_depleted,cr_carbide):
        """ Construction de 3 joints de grain de phases pures
            
            GB_length = longueur de joint simulé en nm
            space_resolution = en nm

            La composition du joint est définie par la liste pos.
            pos[0] est la phase exposée à la surface
            A la profondeur p x space_resolution :
            joint 0 : pos[p] = cr_nominal, la phase est la teneur nominale en Cr
            joint 1 : pos[p] = cr_depleted, la phase est une zone déchromée
            joint 2 : pos[p] = cr_carbide, la phase est un carbure de Cr
        """

        jdg_phases=dict()

        # Création d'un joint de grain tel que la teneur en Cr = cr_nominal
        pos = [cr_nominal]*int(self.GBlength)
        pos2 = []
        for p in range(len(pos)):
            pos2.append([pos[p],pos[p]])

        jdg_phases[0]=pos2

        # Création d'un joint de grain tel que la teneur en Cr = cr_depleted
        pos = [cr_depleted]*int(self.GBlength)
        pos2 = []
        for p in range(len(pos)):
            pos2.append([pos[p],pos[p]])

        jdg_phases[1]=pos2

        # Création d'un joint de grain tel que la teneur en Cr = cr_carbide
        pos = [cr_carbide]*int(self.GBlength)
        pos2 = []
        for p in range(len(pos)):
            pos2.append([pos[p],pos[p]])

        jdg_phases[2]=pos2

        return jdg_phases

    
    def randomGB(self,cr_mean,cr_sd,cr_depleted,cr_carbide):
        """ Tirage aléatoire de la position des carbures et construction des
            phases d'un joint de grain
            
            GB_length = longueur de joint simulé en nm
            space_resolution = en nm
            GBC = taux de couverture du joint de grain par les carbures de Cr
            rC = longueur d'un carbure de Cr

            La composition du joint est définie par la liste pos.
            pos[0] est la phase exposée à la surface
            A la profondeur p x space_resolution :
            pos[p] = cr_nominal si la phase est la teneur nominale en Cr
            pos[p] = cr_depleted si la phase est une zone déchromée
            pos[p] = cr_carbide si la phase est un carbure de Cr

            Un jdg est composé de deux compositions chimiques, l'une initiale,
            l'autre courante (modifiée par la RIS par exemple)
        """
        #print 'nb jdg=',self.nb_gbs
        jdg_phases=dict()

        for my_jdg in range(self.nb_gbs):
            
            # Nombre de carbures à disposer sur chaque joint de grain de longueur GB_length
            nb_carbides_f = self.GBlength*self.GBC/2/self.rC
	    
            # Création d'un joint de grain de teneur nominale en Cr
            pos = int(self.GBlength) * [cr_mean]
               
	    # Random distribution of carbides along grain boundary
            if nb_carbides_f == 0:
                nb_carbides = 0

            else:
                nb_carbides=int(round(nb_carbides_f))
                if nb_carbides_f < 1:
                    if random()<self.GBC:
                        nb_carbides = 1
                    else:
                        nb_carbides = 0
		
                if nb_carbides > 0:
                    # Création de la liste des indices possibles des germes
                    # de carbures
                    indices_germes=[]
                    for ig in range(int(self.rC-1),int(self.GBlength-self.rC+1)):
                        indices_germes.append(ig)

                    # Création des carbures
                    carb = 1
                    while carb <= nb_carbides:

                        # Tirage aléatoire d'un germe dans la liste
                        my_germe = sample(indices_germes,1)[0]

                        # Vérification que l'espace requis pour la croissance
                        # du germe est disponible

                        indices_requis = []
                        for ir in range(my_germe-self.rC,my_germe+self.rC):
                            indices_requis.append(ir)
                        growth = True
                        for i in indices_requis:
                            if (i<0 or i>len(pos)-1):
                                growth = False
                            if (i>=0 and i<=len(pos)-1):
                                if pos[i]==cr_carbide:
                                    growth = False
                                    if i in indices_germes:
                                        indices_germes.remove(i)

                        # Croissance du germe (si possible) et
                        # réduction de la liste des indices de germes
                        if growth:
                            for i in indices_requis:
                                pos[i]=cr_carbide
                        
                        # Depletion en chrome (si zCr>0)
                        if growth:
                            indices_depletion = []
                            for idep in range(my_germe-self.rC-self.zCr,my_germe+self.rC+self.zCr):
                                indices_depletion.append(idep)
                            for i in indices_depletion:
                                if (i>=0 and i<=len(pos)-1):
                                    if pos[i]!=cr_carbide:
                                        pos[i]=cr_depleted

                        # Réduction de la liste des indices de germes
                        if growth:
                            for i in indices_requis:
                                if i in indices_germes:
                                    indices_germes.remove(i)

                        # Incrémentation du carbure
                        if growth:
                            carb+=1
            # duplication du joint
            pos2 = []
            for p in range(len(pos)):
                pos2.append([pos[p],pos[p]])
                            
            jdg_phases[my_jdg]=pos2

        return jdg_phases

    def Coef_T(self):
        """ Calcul de l'effet de la température sur la cinétique d'oxydation.
            material_type : type du matériau mater
            temperature : en K
            oxA : coefficient valable pour le type de matériau
            oxQ : energie d'activation de l'oxydation du matériau
            Coef_T : sans unité """

        R = 8.315
        temperatureK = self.T + 273.15
        coef_temp = exp(-self.oxQ/R/temperatureK)

        return coef_temp 

    def compute_DEcP(self,T,H2):
        """ Calcul de l'écart DEcP au potentiel d'équilibre de Ni/NiO à partir
            de la température et de la teneur en hydrogène dissous dans l'eau.
            DEcP : en mV
            H2 : en ml/kg eau
            T : en C """

        H2_eq_Ni_NiO = 2.e-6 * exp(0.0256*(T+273.15))
        decp = 1000*((8.314*(T+273.15))/(2*96500))*log(H2/H2_eq_Ni_NiO)

        return decp

    def Coef_H(self,decp):
        """ Calcul de l'effet de la teneur en H dissous sur la cinétique d'oxydation.
            decp : écart au potentiel d'équilibre Ni/NiO en mV
            coef_H : sans unité
        """
        coef_H = self.oxg1 + self.oxg2*exp(-self.oxg3*decp)
        
        return coef_H

    def OxidationTime(self, p, XCr,Coef_H,Coef_T,dp):
        """ Temps (h) time_ox necessaire pour un incrément d'oxydation
            de 1 nm.
            time_ox : incrément de temps (h)
            p : profondeur (nm) dejà oxydée
            XCr : teneur (%) en Cr
            ox : matrice des coefficients d'oxydation des phases
            """
        
        # Chargement des coefficients d'oxydation des 3 phases
        a = self.a
        b0 = self.b0 
        b1 = self.b1 
        c0 = self.c0 
        c1 = self.c1 
        
        # construction de la matrice des coefficients d'oxydation
        b = b0 * exp(-b1 * XCr)
        c = c0 * exp(c1 * XCr)
        ox = [a, b * Coef_H * Coef_T, c]

        # Calcul du temps d'oxydation de dp (nm) de la phase en avant du front
        # d'oxydation situé à la profondeur pre_ox

        expi = exp((p - ox[0])/ox[1])
        expf = exp((p + dp - ox[0])/ox[1])
        
        time_ox = (expf - expi) / ox[2]

        return time_ox


    def OxideGrowth(self,pos):
        """ Calcul la cinétique d'oxydation du joint de grain de
            longueur GB_length.

            oxide_growth : liste des couples (temps,profondeur) d'oxydation
            myelectrochem : classe electrochem
            mater : matériau où se trouve le jdg
            t : temps (en h) au bout duquel l'oxyde a atteint la profondeur p
            p : profondeur (en nm) atteinte au temps t
            dp : taille des incréments d'oxydation imposé à 1 nm
            myGB : table des phases (0, 1 ou 2) du jdg où l'oxyde progresse
            pox_c : longueur de joint à oxyder (en nm) pour permettre
                        l'amorçage.
            GB_length : longueur max de joint à oxyder (en nm) pour permettre
                        l'amorçage.
            max_time : temps limite du calcul de CSC

            """
        t = 0
        p = 0
        oxidation_depth = [p]
        oxidation_time = [t]
        dp = 1
        pox_c = int(self.Poxc)
        GB_length = int(self.GBlength)
        max_time = self.stop_time
        decp = self.compute_DEcP(self.T,self.DH)
        coef_H = self.Coef_H(decp)
        coef_T = self.Coef_T()

        while (t <= max_time and p < GB_length):
            # Mise à jour teneur en Cr sous l'effet du RIS
            fluence = self.fluence
            a = self.r1
            b = self.r2
            c = self.r3
            Cr_0 = self.Cr_0

            pos = self.RIS(a,b,c,pos,fluence,Cr_0)
            
            # Calcul du temps nécessaire pour oxyder de p+dp
            t += self.OxidationTime(p,pos[p][1],coef_H,coef_T,dp)
            p += dp
            if not isnan(t):
                oxidation_time.append(t)
                oxidation_depth.append(p)

        return oxidation_time, oxidation_depth

if __name__ == "__main__":

    run_modele = LocalInitiationLaw()

    # Inputs

    language = 'english'
    bins = 100

    # Grain boundaries
    # GBlength, rC, zCr en nm
    # cr en %
    # 0 < GBC < 1
    
    run_modele.GBlength = 1000
    run_modele.nb_gbs = 10
    run_modele.cr_nominal = 16.
    run_modele.cr_sd = 2.
    run_modele.cr_depleted = 12.
    run_modele.cr_carbide = 57.
    
    run_modele.GBC = 0.2
    run_modele.rC = 50
    run_modele.zCr = 10

    # Critical depth to initiation
    run_modele.Poxc = 300

    # Effect of EcP
    run_modele.oxg1 = 0.1
    run_modele.oxg2 = 0.3
    run_modele.oxg3 = 0.04

    # Effect of temperature
    run_modele.oxQ = 139000.

    # Cr effect
    run_modele.a = 0.
    run_modele.b0 = 8.67e14
    run_modele.b1 = 6.22e-2
    run_modele.c0 = 3.e-3
    run_modele.c1 = 7.5e-2

    # RIS
    run_modele.r1 = -0.5
    run_modele.r2 = 0.5
    run_modele.r3 = 0.2
    run_modele.Cr_0 = 16.

    # Environmental parameters
    run_modele.T = 360.
    run_modele.DH = 30.
    run_modele.fluence = 0.

    # Stopping criterion
    run_modele.stop_time = 1.e4

    # Time to check the frequency of oxidation depths
    # time_frequency < stop_time
    time_frequency = 1000. 

    if time_frequency > run_modele.stop_time:
        time_frequency = run_modele.stop_time
        print("Warning: selected time_frequency is larger than stop_time")
              
    # Création des jdg
    pos = run_modele.randomGB(run_modele.cr_nominal,run_modele.cr_sd,run_modele.cr_depleted,run_modele.cr_carbide)
    
    # Création de 3 jdg de phases pures différentes
    pure_pos = run_modele.pureGB(run_modele.cr_nominal,run_modele.cr_depleted,run_modele.cr_carbide)
    
    # Impression des jdg
    run_modele.PrintGbs(pos,run_modele.cr_nominal,run_modele.cr_depleted,run_modele.cr_carbide)

    # Calcul de la cinétique d'oxydation des phases pures de jdg, sans effet
    # du RIS. Si une fluence est imposée, elle est temporairement annulée.
    pure_oxide_growth = []
    my_fluence = run_modele.fluence
    run_modele.fluence = 0.
    for j in range(len(pure_pos)):
        pure_oxide_growth.append(run_modele.OxideGrowth(pure_pos[j]))

    run_modele.PlotPureOxidationKinetics(pure_oxide_growth,language)

    # Trace l'élolution de la teneur en Cr des phases pures en fonction de la fluence
    run_modele.fluence = my_fluence
    run_modele.PlotRIS('english') 
    
    # Calcul de la cinétique d'oxydation des jdg
    oxide_growth = []

    for j in range(len(pos)):
        oxide_growth.append(run_modele.OxideGrowth(pos[j]))

    run_modele.PlotOxidationKinetics(pure_oxide_growth,oxide_growth,language)

