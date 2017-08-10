import algo
import matplotlib.pyplot as pl

###Choose what you want to do : 0 = no, 1 = yes ###
execute = [True,True] # [given configuration, real time]
detail = True      # 1 : markers added step by step (only used in execute[0])

### Parameters ###
a = 1 #width
b = 1 #height
h_ini = 1/4 #1/2, 1/4 or 1/6 for example         # the lower it is the more we have markers and cells at the beginning
delta_t = 1.E-2 #1.E-2   # time step
nb_iter = 1000  #1000    # total number of steps

### Variables ###
h = h_ini/2
generation_old = 0
generation = 0
compteur = 0             # current step number
dist_min = 0

### Markers ###
marker1 = 'go'           # green circle markers : initial positions of added markers
marker2 = 'yo'           # green yellow markers : initial positions of initial markers
marker3 = 'r+'           # red + markers : current positions of added markers
marker4 = 'k+'           # black + markers : current position of initial markers

### Initialization ###
stock_gen = [0]
M_ini = [algo.initialisation(a,b,h_ini,0,generation,stock_gen,1)]
stock_gen = [0]
M = algo.initialisation(a,b,h_ini,0,generation,stock_gen,1)
C = algo.initialisation(a,b,h_ini,1,generation,None,0)
generation = stock_gen[0]

C_n = C.shape[0]
C_m = C.shape[1]
M_n = M_ini[0].shape[0]
M_m = M_ini[0].shape[1]
nb_ini = M_n*M_m
M_ini_utile = []
M_utile = []
for m_j in range(M_m):
    for m_i in range(M_n):
        marqueur_new_ini = algo.marqueur(0,0,0,0,-1,0)
        marqueur_new_ini.copie(M_ini[0][m_i,m_j])
        M_ini_utile.append(marqueur_new_ini)
        M_utile.append(M[m_i,m_j])

n = 2*int(a/h_ini)+1
m = 2*int(b/h_ini)+1
M_ini2 = algo.initialisation(a,b,h_ini/2,0,-1,None,0)
for j in range(m):
    for i in range(n):
        if M_ini[-1][i,j].generation >= 0:
            M_ini2[2*i,2*j].generation = M_ini[-1][i,j].generation
M_ini.append(M_ini2)
for e in range(len(M_utile)):
    M_ini_utile[e].i = 2*M_ini_utile[e].i
    M_ini_utile[e].j = 2*M_ini_utile[e].j
    M_utile[e].i = 2*M_utile[e].i
    M_utile[e].j = 2*M_utile[e].j
variables = [M_ini,h,generation]
generation_old = generation

### Given configuration ###
### We do nb_iter moves and then we apply the algorithm to the configuration we got ###
if execute[0]:
    compteur = nb_iter
    ### nb_iter moves ###
    algo.mouvement_liste(M_utile,compteur,delta_t)
    tri = algo.tri_liste(M_utile,M_ini_utile,nb_ini)
    maillage_ini = algo.maillage_ini(C)
    algo.affiche(tri,maillage_ini,C,h_ini,1,None,0)
    ### add markers ###
    for c_i in range(C_n):
        for c_j in range(C_m):
            test_distance = False
            M_ini_hd = algo.marqueur(0,0,0,0,-1,0)
            M_ini_bd = algo.marqueur(0,0,0,0,-1,0)
            M_ini_bg = algo.marqueur(0, 0, 0, 0, -1, 0)
            M_ini_hg = algo.marqueur(0, 0, 0, 0, -1, 0)
            dist_hd = 1000
            dist_bd = 1000
            dist_bg = 1000
            dist_hg = 1000
            for m in range(len(M_utile)):
                if C[c_i,c_j].distance(M_utile[m]) <= h_ini/2:
                    test_distance = True
                    break
                else:
                    if M_utile[m].x >= C[c_i,c_j].x and M_utile[m].y >= C[c_i,c_j].y:
                        distance = algo.distance_cellule(C[c_i,c_j].x,C[c_i,c_j].y,M_utile[m].x,M_utile[m].y,h_ini,0)
                        if distance < dist_hd and distance > h_ini/8:
                            M_ini_hd.copie(M_ini_utile[m])
                            dist_hd = distance
                    elif M_utile[m].x >= C[c_i,c_j].x and M_utile[m].y <= C[c_i,c_j].y:
                        distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x, M_utile[m].y,
                                                         h_ini, 1)
                        if distance < dist_bd and distance > h_ini/8:
                            M_ini_bd.copie(M_ini_utile[m])
                            dist_bd = distance
                    elif M_utile[m].x <= C[c_i, c_j].x and M_utile[m].y <= C[c_i, c_j].y:
                        distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x, M_utile[m].y,
                                                         h_ini, 2)
                        if distance < dist_bg and distance > h_ini/8:
                            M_ini_bg.copie(M_ini_utile[m])
                            dist_bg = distance
                    elif M_utile[m].x <= C[c_i, c_j].x and M_utile[m].y >= C[c_i, c_j].y:
                        distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x, M_utile[m].y,
                                                         h_ini, 3)
                        if distance < dist_hg and distance > h_ini/8:
                            M_ini_hg.copie(M_ini_utile[m])
                            dist_hg = distance
            if test_distance == False: ### if no marker at distance within h_ini/2 : ###
                variables = algo.remaillage(C[c_i,c_j],h,h_ini,delta_t,compteur,M_ini,M_ini_utile,M_utile,a,b,generation)
                M_ini = variables[0]
                h = variables[1]
                generation = variables[2]
                marqueur_test = variables[3]
                ### add neighbours for finite differences ###
                if detail: # if detail = True : add marker one by one
                                # The pink cross is the exact position computed for the new marker
                    tri = algo.tri_liste(M_utile,M_ini_utile,nb_ini)
                    algo.affiche(tri,maillage_ini,C,h_ini,1,marqueur_test,1)
                generation_old = generation
    ### Now we display the result ###
    tri = algo.tri_liste(M_utile,M_ini_utile,nb_ini)
    algo.affiche(tri,maillage_ini,C,h_ini,1,None,0)
    
### real time ###
if execute[1]:
    fig, ax = pl.subplots()
    tri = algo.tri_liste(M_utile,M_ini_utile,nb_ini)
    algo.tracecarre(C,h_ini)
    while compteur <= nb_iter: # we want nb_iter steps
        if compteur == 0:
            maillage_ini = algo.maillage_ini(C)
            points, = ax.plot(tri[4], tri[5], marker1, linestyle='None')
            points2, = ax.plot(maillage_ini[0], maillage_ini[1], marker2, linestyle='None')
            points3, = ax.plot(tri[2], tri[3], marker3, linestyle='None')
            points4, = ax.plot(tri[0], tri[1], marker4, linestyle='None')
            ax.set_xlim(-1.5, 1.5) 
            ax.set_ylim(-1.5, 1.5)
        else:
            algo.mouvement_liste(M_utile,1,delta_t)
            new_points = False
            for c_i in range(C_n):
                for c_j in range(C_m):
                    test_distance = False
                    M_ini_hd = algo.marqueur(0, 0, 0, 0, -1, 0)
                    M_ini_bd = algo.marqueur(0, 0, 0, 0, -1, 0)
                    M_ini_bg = algo.marqueur(0, 0, 0, 0, -1, 0)
                    M_ini_hg = algo.marqueur(0, 0, 0, 0, -1, 0)
                    dist_hd = 1000
                    dist_bd = 1000
                    dist_bg = 1000
                    dist_hg = 1000
                    for m in range(len(M_utile)):
                        if C[c_i, c_j].distance(M_utile[m]) <= h_ini / 2:
                            test_distance = True
                            break
                        else:
                            if M_utile[m].x >= C[c_i, c_j].x and M_utile[m].y >= C[c_i, c_j].y:
                                distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x,
                                                                 M_utile[m].y, h_ini, 0)
                                if distance < dist_hd and distance > dist_min:
                                    M_ini_hd.copie(M_ini_utile[m])
                                    dist_hd = distance
                            elif M_utile[m].x >= C[c_i, c_j].x and M_utile[m].y <= C[c_i, c_j].y:
                                distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x,
                                                                 M_utile[m].y,
                                                                 h_ini, 1)
                                if distance < dist_bd and distance > dist_min:
                                    M_ini_bd.copie(M_ini_utile[m])
                                    dist_bd = distance
                            elif M_utile[m].x <= C[c_i, c_j].x and M_utile[m].y <= C[c_i, c_j].y:
                                distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x,
                                                                 M_utile[m].y,
                                                                 h_ini, 2)
                                if distance < dist_bg and distance > dist_min:
                                    M_ini_bg.copie(M_ini_utile[m])
                                    dist_bg = distance
                            elif M_utile[m].x <= C[c_i, c_j].x and M_utile[m].y >= C[c_i, c_j].y:
                                distance = algo.distance_cellule(C[c_i, c_j].x, C[c_i, c_j].y, M_utile[m].x,
                                                                 M_utile[m].y,
                                                                 h_ini, 3)
                                if distance < dist_hg and distance > dist_min:
                                    M_ini_hg.copie(M_ini_utile[m])
                                    dist_hg = distance
                    if test_distance == False: ### if no marker within h_ini/2 : ###
                        new_points = True
                        variables = algo.remaillage2(C[c_i, c_j], h, h_ini, delta_t, compteur, M_ini, M_ini_utile,
                                                     M_utile, a, b,
                                                     generation, M_ini_hd, M_ini_bd, M_ini_bg, M_ini_hg, dist_hd,
                                                     dist_bd, dist_bg, dist_hg)
                        M_ini = variables[0]
                        h = variables[1]
                        generation = variables[2]
                        marqueur_test = variables[3]
                        test = algo.test_active_voisin(M_ini_utile,h_ini,generation_old)
                        if test[0] == False:
                            algo.active_voisin_horizontal(marqueur_test,M_ini_utile[generation_old],M_ini,M_ini_utile,M_utile,a,b,h_ini,h,compteur,delta_t,variables)
                            generation = variables[2]
                        if test[1] == False:
                            algo.active_voisin_vertical(marqueur_test,M_ini_utile[generation_old],M_ini,M_ini_utile,M_utile,a,b,h_ini,h,compteur,delta_t,variables)
                            generation = variables[2]
                        generation_old = generation
            tri = algo.tri_liste(M_utile,M_ini_utile,nb_ini)
            if new_points == True:
                points.set_data(tri[4], tri[5])
            points3.set_data(tri[2], tri[3])
            points4.set_data(tri[0], tri[1])
        pl.pause(0.001)
        compteur+=1
    pl.show()
    if compteur == nb_iter:
        print("fini")