"""
Les modules suivants doivent etre installes:
sys
pandas
numpy
os
scipy
matplotlib

Pour installer: 
-'pip install nom_module'

Importation des fonctions de traitement des donnees
help(nom_du_sous_module.nom_fonction) donne la notice de la fonction
Les fonction sont acessibles et modifiables dans 
exploitation_signaux/nom_du_sous_module

"""
import exploitation_signaux.extraction_data_set as ex
import exploitation_signaux.traitement_data_set as td
import exploitation_signaux.analyse_data_set as ad

from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np

    

"""
Recuperation des signaux d'interet avec la fonction d'extraction situee dans 
exploitation_signaux.extraction_data_set il faut fournir en premier le chemin
du dossier ou sont stockes les signaux et en deuxieme le chemin vers le fichier
de calibration
"""
H,Im1,Re1,Liste_Frequences,Courant=ex.extraction(r"D:\IMT-Atlantique\1Année\CODEV\DATA\Mesure_FMR_Co60nm - Copie",r"D:\IMT-Atlantique\1Année\CODEV\DATA\Etalonnage electro-magnet\Etalonage bis (2).txt")


def Lorentzienne(f,A,Df,fres):
    return A*Df**2/((f-fres)**2+Df**2)

def AsLorentzienne(f,A,Df,fres):
    return A*(Df)*(f-fres)/((f-fres)**2+Df**2)

def ClDeLorentzienne(f,amplitude,Df,fres,a,b):
    return a*Lorentzienne(f,amplitude,Df,fres_empirique)+b*AsLorentzienne(f,amplitude,Df,fres_empirique)

"""
Listes de stockage des parametres de la combinaison lineaire a*Lorentzienne+b*AsLorentzienne
"""
A=[]
B=[]


"""
La premiere liste contiendra en position i l'amplitude du signal Im[i],
 la deuxieme la demi largeur a mi hauteur, 
 la 3 eme la frequence de resonnance
 la 4 eme la frequence estimee du second pic
 """
Listes_Des_Parametre_De_Fit=[[],[],[],[]] 


n=len(Im1)

"""Definition de la largeur (en nombre de points de la taille de la fenetre)"""
taille_fenetre=50



for i in range(n):
    
    nparray_Im_De_Travail = Im1[i]
    Liste_Im_Fenetree, Xmodf = td.fenetrage(nparray_Im_De_Travail[:630], Liste_Frequences[:630], taille_fenetre)
    Liste_Im_De_Travail = list(nparray_Im_De_Travail)

    """Calcul du second pic"""
    ind_fres = Liste_Im_De_Travail.index(max(Liste_Im_De_Travail[:630]))
    fres_empirique = Liste_Frequences[ind_fres]
    v, ind_fres_sec = ad.secondPic(Liste_Im_De_Travail, 683, ind_fres, 20, 9)
    if v:
        Listes_Des_Parametre_De_Fit[3].append(Liste_Frequences[ind_fres_sec])
        plt.scatter(Liste_Frequences[ind_fres_sec],
        nparray_Im_De_Travail[ind_fres_sec], color='red')
    else:
        Listes_Des_Parametre_De_Fit[3].append(-1)

    """Affichage dans le terminal de l'avancement"""
    td.update_progress(int(i*100/(n-1)))
    # print("=> "+str(int(i*100/(n-1)))+"%",end='\r')

    """
    Premiere estimation des parametres amplitude ,demi largeur a mi hauteur,
    frequence de resonnance
    """
    params = [max(Liste_Im_Fenetree), 1, fres_empirique]
    if Listes_Des_Parametre_De_Fit[1] != []:
        params[1] = Listes_Des_Parametre_De_Fit[1][-1]
    amplitude = params[0]

    """
    
    """
    parametres2, cov = curve_fit(ClDeLorentzienne, Xmodf, Liste_Im_Fenetree, p0=[params[0], params[1], params[2],  0.707, 0.707])

    Listes_Des_Parametre_De_Fit[0].append(parametres2[0])
    Listes_Des_Parametre_De_Fit[1].append(parametres2[1])
    Listes_Des_Parametre_De_Fit[2].append(parametres2[2])
    A.append(parametres2[3])
    B.append(parametres2[4])

    """Trace et sauvegarde des figures"""
    plt.plot(Liste_Frequences, nparray_Im_De_Travail, color='blue')
    plt.plot(Liste_Frequences, ClDeLorentzienne(
    Liste_Frequences, parametres2[0], parametres2[1], parametres2[2], parametres2[3], parametres2[4]), color='red')
    plt.title(str(H[i])+" Tesla ; Courant de : "+str(Courant[i])+"A")
    """Répertoire à spécifier : """
    plt.savefig(r"D:\IMT-Atlantique\1Année\CODEV\Courbes\Partie_Imaginaire_CL2\\" +str(i)+"_"+str(H[i])+"Tesla;Courant_de_"+str(Courant[i])+"A"+".png")
    plt.clf()



"""Trace des coefficients des coefficients trouves par le fit de la CL de l'antiLorentzienn
et de le Lorentzienne"""
# plt.plot(np.array(A)**2+np.array(B)**2)



"""Trace des parametres du fit en fonction de H"""

"""fres(n=1) en fonction de H"""
plt.figure(1)
Frequence_Second_pic=Listes_Des_Parametre_De_Fit[3]
H_avec_second_pic_detecte=[]
FresPSSW=[]
for i in range(len(Frequence_Second_pic)):
    if Frequence_Second_pic[i]!=-1:
        H_avec_second_pic_detecte.append(H[i])
        FresPSSW.append(Frequence_Second_pic[i])
plt.plot(H_avec_second_pic_detecte,FresPSSW)
plt.title("fres(n=1) (FresPSSW) en fonction de H")
plt.show()

"""fres(n=0) en fonction de H"""
plt.figure(2)
Frequences_Resonnance_mode_fondamental=Listes_Des_Parametre_De_Fit[2]
plt.plot(H,Frequences_Resonnance_mode_fondamental)
plt.title("fres(n=0) (FresFMR) en fonction de H")
plt.show()

"""Amplitude en fonction de H"""
plt.figure(3)
Amplitude=Listes_Des_Parametre_De_Fit[0]
plt.plot(H,Amplitude)
plt.title("Amplitude en fonction de H")
plt.show()


"""Code sur le Df"""
plt.figure(4)
plt.plot(H,Listes_Des_Parametre_De_Fit[1],color='blue')
# plt.plot(Listes_Des_Parametre_De_Fit[2], Listes_Des_Parametre_De_Fit[1], 'b-', label='data', color='green')
for i in range(len(Listes_Des_Parametre_De_Fit[1])):
    plt.scatter(H[i],Listes_Des_Parametre_De_Fit[1][i],color='red')

# slope1, intercept1, r_value1, p_value1, std_err1 = linregress(Hprim[:247], LargeurMiHauteur[:247])
# slope2, intercept2, r_value2, p_value2, std_err2 = linregress(Hprim[247:], LargeurMiHauteur[247:])
# plt.plot(Hprim[247:],np.array(Hprim[247:])*slope2+intercept2,color='blue')
# plt.plot(Hprim[:247],np.array(Hprim[:247])*slope1+intercept1,color='blue')


plt.title("Demi Largeur a mi hauteur en focntion du champ exterieur H")
plt.show()
 
"""Partie sur la difference des carres des frequences de 
resonnance du premier mode et du second mode"""  

plt.figure(5)
plt.title("f²(n=1)-f²(n=0) en fonction de H")
H_avec_second_pic_detecte = []
deltaFresCarre = []
Frequence_Second_pic = Listes_Des_Parametre_De_Fit[3]
Frequences_Resonnance_mode_fondamental = Listes_Des_Parametre_De_Fit[2]

for i in range(len(Frequence_Second_pic)):
    """Si il y a un second pic alors on recupere le champ
    et on calcule la difference des carres"""
    if Frequence_Second_pic[i] != -1:
        
        H_avec_second_pic_detecte.append(H[i])
        deltaFresCarre.append((Frequence_Second_pic[i]**2-Frequences_Resonnance_mode_fondamental[i]**2))
        plt.scatter(H[i], Frequence_Second_pic[i]**2 - Frequences_Resonnance_mode_fondamental[i]**2, color='red')

H_avec_second_pic_detecte_abs = [abs(x) for x in H_avec_second_pic_detecte]

indice_coupure = H_avec_second_pic_detecte_abs.index(min(H_avec_second_pic_detecte_abs))

slope1, intercept1, r_value1, p_value1, std_err1 = linregress( H_avec_second_pic_detecte[:indice_coupure], deltaFresCarre[:indice_coupure])
slope2, intercept2, r_value2, p_value2, std_err2 = linregress( H_avec_second_pic_detecte[indice_coupure:], deltaFresCarre[indice_coupure:])

plt.plot(H_avec_second_pic_detecte[indice_coupure:], np.array( H_avec_second_pic_detecte[indice_coupure:])*slope2+intercept2, color='blue')
plt.plot(H_avec_second_pic_detecte[:indice_coupure], np.array( H_avec_second_pic_detecte[:indice_coupure])*slope1+intercept1, color='blue')


plt.show()



"""Fit pour la formule de Kittel"""
mu0=12.566370*10**(-7)
def kittel(x, gamma, meff):
    return (gamma)*(x*(x+meff))**0.5
Hfresabs=[abs(x) for x in H]

"""Version champ positif"""
plt.figure(6)

indice_coupure_bas=500
indice_coupure_haut=Hfresabs.index(min(Hfresabs))+4

paramsk, covk = curve_fit(kittel, np.array(H[indice_coupure_haut:indice_coupure_bas]), Frequences_Resonnance_mode_fondamental[indice_coupure_haut:indice_coupure_bas], p0=[30000000 , 1.7])
plt.plot(np.array(H[indice_coupure_haut:indice_coupure_bas]), Frequences_Resonnance_mode_fondamental[indice_coupure_haut:indice_coupure_bas], 'b-', label='data')
plt.plot(np.array(H[indice_coupure_haut:indice_coupure_bas]), kittel(np.array(H[indice_coupure_haut:indice_coupure_bas]), paramsk[0], paramsk[1] ), 'r-', label='fit')
plt.xlabel("Champ extérieur Hext en Tesla")
plt.label("f(n=0)")
plt.title("f(n=0) en fonction de H; gamma = "+str(paramsk[0])+"GHz/Tesla ; Meff = " + str(paramsk[1]))
plt.show()

"""Version champ négatif"""
plt.figure(7)

indice_coupure_bas=Hfresabs.index(min(Hfresabs))-4
indice_coupure_haut=0
paramsk, covk = curve_fit(kittel, -np.array(H[indice_coupure_haut:indice_coupure_bas]), Frequences_Resonnance_mode_fondamental[indice_coupure_haut:indice_coupure_bas], p0=[36 , 1.7])
plt.plot(-np.array(H[indice_coupure_haut:indice_coupure_bas]), Frequences_Resonnance_mode_fondamental[indice_coupure_haut:indice_coupure_bas], 'b-', label='data')
plt.plot(-np.array(H[indice_coupure_haut:indice_coupure_bas]), kittel(-np.array(H[indice_coupure_haut:indice_coupure_bas]), paramsk[0], paramsk[1] ), 'r-', label='fit')
plt.title("f(n=0) en fonction de H; gamma = "+str(paramsk[0])+"GHz/Tesla ; Meff = " + str(paramsk[1]))
plt.show()


                  