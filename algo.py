import numpy as np
import matplotlib.pyplot as plt

### marqueur class ###
class marqueur:
    def __init__(self,x,y,i,j,generation,h):
        self.x = x*h
        self.y = y*h
        self.i = i
        self.j = j
        self.generation = generation

    def copie(self,M):
        self.x = M.x
        self.y = M.y
        self.i = M.i
        self.j = M.j
        self.generation = M.generation

    def norme(self):
        return (self.x*self.x + self.y*self.y)

    def distance(self,M):
        X = np.abs(self.x-M.x)
        Y = np.abs(self.y-M.y)
        dist = max(X,Y)
        return dist

    ### return the closest neighbour among M1,M2,M3,M4 ###
    def voisins(self,M1,M2,M3,M4):
        proche = [M1,M2,M3,M4]
        for e in range(3):
            for f in range(4):
                if e < f and self.distance(proche[e]) > self.distance(proche[f]):
                    temp = proche[e]
                    proche[e] = proche[f]
                    proche[f] = temp
        return proche

    def vitesse(self): # Velocity field
        if self.norme() > 1:
            vit = [0,0]
        else:
            vit = [(-1+np.sqrt(self.norme()))*self.y,-(-1+np.sqrt(self.norme()))*self.x]
        return vit

    def deplacement(self,delta_t): # Move forward
        vit = self.vitesse()
        self.x += vit[0]*delta_t
        self.y += vit[1]*delta_t

    def deplacement_arriere(self,delta_t): # Move backward
        vit = self.vitesse()
        self.x -= vit[0]*delta_t
        self.y -= vit[1]*delta_t

### cellule class ###
class cellule:
    def __init__(self,i,j,h):
        self.x = i*h
        self.y = j*h
        self.i = i
        self.j = j
        self.h = h

    def distance(self,M):
        X = np.abs(self.x-M.x)
        Y = np.abs(self.y-M.y)
        dist = max(X,Y)
        return dist


def maillage_ini(C):
    C_x = []
    C_y = []
    q = C.shape[0]
    p = C.shape[1]
    for i in range(p):
        for j in range(q):
            C_x.append(C[i, j].x)
            C_y.append(C[i, j].y)
    return [C_x, C_y]

### Initialization of markers array###
def initialisation(a,b,h,classe,generation,stock_gen,mode):
    n = 2*int(a/h)+1
    m = 2*int(b/h)+1
    if classe == 0:
        M = np.eye(n,m,dtype=marqueur)
    else:
        M = np.eye(n,m,dtype=cellule)
    for j in range(m):
        v = int((2*j-m+1)/2)
        for i in range(n):
            u = int((n-1-2*i)/2)
            if classe == 0:
                M[i,j] = marqueur(v,u,i,j,generation,h)
                if mode == 1:
                    generation+=1
            else:
                M[i,j] = cellule(v,u,h)
    if classe ==0 and mode == 1:
        stock_gen[0] = generation
    return M

### Display ###
def tracecarre(C,h):
    n = C.shape[0]
    m = C.shape[1]
    for i in range(n):
        for j in range(m):
            x = []
            y = []
            
            x.append(C[i,j].x-h/2)
            y.append(C[i,j].y-h/2)
            
            x.append(C[i,j].x+h/2)
            y.append(C[i,j].y-h/2)

            x.append(C[i,j].x+h/2)
            y.append(C[i,j].y+h/2)

            x.append(C[i,j].x-h/2)
            y.append(C[i,j].y+h/2)

            x.append(C[i,j].x-h/2)
            y.append(C[i,j].y-h/2)

            plt.plot(x,y,'b-')

def affiche(tri,maillage_ini,C,h,mode,marqueur_test,mode2):
    plt.plot(tri[4], tri[5], 'go')
    plt.plot(maillage_ini[0], maillage_ini[1], 'yo')
    if mode == 1:
        plt.plot(tri[2], tri[3], 'r+')
        plt.plot(tri[0], tri[1], 'k+')
    if mode2 == 1:
        plt.plot(marqueur_test.x,marqueur_test.y, 'mx')
    tracecarre(C,h)
    plt.axis([-1.2,1.2,-1.2,1.2])
    plt.show()

### Sort markers ###
def tri_liste(M,M_ini,nb_ini):
    M0_x = []
    M0_y = []
    M_newx = []
    M_newy = []
    M_ini_newx = []
    M_ini_newy = []
    for i in range(len(M)):
        if M[i].generation < nb_ini:
            M0_x.append(M[i].x)
            M0_y.append(M[i].y)
        elif M[i].generation >= nb_ini:
            M_newx.append(M[i].x)
            M_newy.append(M[i].y)
            M_ini_newx.append(M_ini[i].x)
            M_ini_newy.append(M_ini[i].y)
    tri = [M0_x,M0_y,M_newx,M_newy,M_ini_newx,M_ini_newy]
    return tri
                
### Moves ###
def mouvement_liste(M,nb_iter,delta_t):
    for k in range(nb_iter):
        for i in range(len(M)):
            M[i].deplacement(delta_t)

### Test if we have markers close enough to compute finite differences ###
def test_active_voisin(M_ini_utile,h_ini,indice):
    horizontal = False
    vertical = False
    for e in range(len(M_ini_utile)):
        if e != indice:
            if horizontal == False:
                if M_ini_utile[e].y == M_ini_utile[indice].y:
                    if abs(M_ini_utile[e].x - M_ini_utile[indice].x) <= h_ini:
                        horizontal = True
            if vertical == False:
                if M_ini_utile[e].x == M_ini_utile[indice].x:
                    if abs(M_ini_utile[e].y - M_ini_utile[indice].y) <= h_ini:
                        vertical = True
            if horizontal == True and vertical == True:
                break
    return [horizontal,vertical]

### Create a right horizontal neighbour with abscissa multiple of h_ini ###
def calcule_j(marqueur_test,marqueur_ini,h_ini,h,b):
    level = int(round(h_ini/h))
    if marqueur_ini.x < 0:
        if marqueur_test.x > int(marqueur_ini.x/h_ini)*h_ini -h_ini/2:
            j = int(marqueur_ini.x/h_ini)*level +int(b/h)
        else:
            j = int(marqueur_ini.x/h_ini -1)*level +int(b/h)
        if np.abs(int(marqueur_ini.x/h_ini) - marqueur_ini.x/h_ini) <= 1.E-5:
            if marqueur_test.x > int(marqueur_ini.x/h_ini)*h_ini -h_ini/2: 
                j = (int(marqueur_ini.x/h_ini)+1)*level +int(b/h)
            else:
                j = (int(marqueur_ini.x/h_ini)-1)*level +int(b/h)
    else:
        if marqueur_test.x > int(marqueur_ini.x/h_ini)*h_ini + h_ini/2:
            j = int(marqueur_ini.x/h_ini +1)*level +int(b/h)
        else:
            j = int(marqueur_ini.x/h_ini)*level +int(b/h)
        if np.abs(int(marqueur_ini.x/h_ini) - marqueur_ini.x/h_ini) <= 1.E-5:
            if marqueur_test.x > int(marqueur_ini.x/h_ini)*h_ini + h_ini/2: 
                j = (int(marqueur_ini.x/h_ini)+1)*level +int(b/h)
            else:
                j = (int(marqueur_ini.x/h_ini)-1)*level +int(b/h)
    if np.abs(j - int(j) -1) <= 1.E-5: ### Avoid bug due to computer roundsT ###
        j+=1
    return j

def active_voisin_horizontal(marqueur_test,marqueur_ini,M_ini,M_ini_utile,M_utile,a,b,h_ini,h,compteur,delta_t,variables):
    i = int(a/h)-int(marqueur_ini.y/h)
    j = calcule_j(marqueur_test,marqueur_ini,h_ini,h,b)
    ### Tests pour éviter les bugs dus aux arrondis machine ###
    if np.abs(M_ini[-1][i,j].y-marqueur_ini.y-h) < 1.E-5:
        i+=1
    elif np.abs(M_ini[-1][i,j].y-marqueur_ini.y+h) < 1.E-5:
        i-=1
    M_ini[-1][i,j].generation = variables[2]
    marqueur_new_ini = marqueur(0,0,0,0,-1,0)
    marqueur_new_ini.copie(M_ini[-1][i,j])
    M_ini_utile.append(marqueur_new_ini)
    marqueur_new = marqueur(0,0,0,0,-1,0)
    marqueur_new.copie(M_ini[-1][i,j])
    for e in range(compteur):
        marqueur_new.deplacement(delta_t)
    M_utile.append(marqueur_new)
    variables[2]+=1

### Create a down vertical neighbour with intercept multiple of h_ini ###
def calcule_i(marqueur_test,marqueur_ini,h_ini,h,a):
    #print(marqueur_ini.y/h_ini)
    level = int(round(h_ini / h))
    if marqueur_ini.y > 0:
        if marqueur_test.y < int(marqueur_ini.y/h_ini)*h_ini +h_ini/2:
            i = -int(marqueur_ini.y/h_ini)*level + int(a/h)
        else:
            i = -int(marqueur_ini.y/h_ini +1)*level + int(a/h)
        if np.abs(int(marqueur_ini.y/h_ini) - marqueur_ini.y/h_ini) <= 1.E-5:
            if marqueur_test.y < int(marqueur_ini.y/h_ini)*h_ini +h_ini/2:
                i = -(int(marqueur_ini.y/h_ini)-1)*level + int(a/h)
            else:
                i = -(int(marqueur_ini.y/h_ini)+1)*level + int(a/h)
    else:
        if marqueur_test.y < int(marqueur_ini.y/h_ini)*h_ini -h_ini/2:
            i = -int(marqueur_ini.y/h_ini -1)*level + int(a/h)
        else:
            i = -int(marqueur_ini.y/h_ini)*level + int(a/h)
        if np.abs(int(marqueur_ini.y/h_ini) - marqueur_ini.y/h_ini) <= 1.E-5:
            if marqueur_test.y < int(marqueur_ini.y/h_ini)*h_ini +h_ini/2:
                i = -(int(marqueur_ini.y/h_ini)-1)*level + int(a/h)
            else:
                i = -(int(marqueur_ini.y/h_ini)+1)*level + int(a/h)
    if np.abs(i - int(i) -1) <= 1.E-14:
        i+=1
    return i

def active_voisin_vertical(marqueur_test,marqueur_ini,M_ini,M_ini_utile,M_utile,a,b,h_ini,h,compteur,delta_t,variables):
    j = int(marqueur_ini.x/h)+int(b/h)
    i = calcule_i(marqueur_test,marqueur_ini,h_ini,h,a)
    if np.abs(M_ini[-1][i,j].x-marqueur_ini.x-h) < 1.E-5:
        j-=1
    elif np.abs(M_ini[-1][i,j].x-marqueur_ini.x+h) < 1.E-5:
        j+=1
    M_ini[-1][i,j].generation = variables[2]
    marqueur_new_ini = marqueur(0,0,0,0,-1,0)
    marqueur_new_ini.copie(M_ini[-1][i,j])
    M_ini_utile.append(marqueur_new_ini)
    marqueur_new = marqueur(0,0,0,0,-1,0)
    marqueur_new.copie(M_ini[-1][i,j])
    for e in range(compteur):
        marqueur_new.deplacement(delta_t)
    M_utile.append(marqueur_new)
    variables[2]+=1

### Get the 4 closest neighbours of marqueur ###
def recup_voisins(M_ini,marqueur,h_ini,h,nb_ini,a,b):
    i= marqueur.i
    j= marqueur.j
    nb = int(round(h_ini/h))
    n = M_ini[-1].shape[0]
    voisins = [None,None,None,None]
    test = [False,False,False,False]
    
    for f in range(nb):
        #up:
        if i-f-1 >= 0 and test[0] == False:
            u = i-f-1
            if M_ini[-1][u,j].generation >=0:
                voisins[0] = M_ini[-1][u,j]
                test[0] = True
        #down:
        if i+f+1 < n and test[1] == False:
            u = i+f+1
            if M_ini[-1][u,j].generation >=0:
                voisins[1] = M_ini[-1][u,j]
                test[1] = True
        #left:
        if j-f-1 >= 0 and test[2] == False:
            v = j-f-1
            if M_ini[-1][i,v].generation >=0:
                voisins[2] = M_ini[-1][i,v]
                test[2] = True
        #right:
        if j+f+1 < n and test[3] == False:
            v = j+f+1
            if M_ini[-1][i,v].generation >=0:
                voisins[3] = M_ini[-1][i,v]
                test[3] = True
        if test == [True,True,True,True]:
            break
    return voisins

### finite difference ###
def difference_finie(marqueur,voisins,M_utile,delta_t,compteur,exact,e):
    differences = [None,None,None,None]
    d0 = np.inf
    d1 = np.inf
    d2 = np.inf
    d3 = np.inf
    g = marqueur.generation

    if exact == 1:
        new_marqueur = flot_forward(marqueur,delta_t,compteur)
    
    if voisins[0] != None:
        g0 = voisins[0].generation
        d0 = marqueur.distance(voisins[0])
        if exact == 1:
            new_voisin_h = flot_forward(voisins[0],delta_t,compteur)
    if voisins[1] != None:
        d1 = marqueur.distance(voisins[1])
        g1 = voisins[1].generation
        if exact == 1:
            new_voisin_b = flot_forward(voisins[1],delta_t,compteur)
    if voisins[2] != None:
        d2 = marqueur.distance(voisins[2])
        g2 = voisins[2].generation
        if exact == 1:
            new_voisin_g = flot_forward(voisins[2],delta_t,compteur)
    if voisins[3] != None:
        d3 = marqueur.distance(voisins[3])
        g3 = voisins[3].generation
        if exact == 1:
            new_voisin_d = flot_forward(voisins[3],delta_t,compteur)

    arrondi = 15
    #vertical
    if d0 == d1:
        if exact == 0:
            differences[0] = (round(M_utile[g1].x,arrondi) -round(M_utile[g0].x,arrondi))/(2*d0)
            differences[1] = (round(M_utile[g1].y,arrondi) - round(M_utile[g0].y,arrondi))/(2*d0)
        else:
            differences[0] = (round(new_voisin_b[0],arrondi) - round(new_voisin_h[0],arrondi))/(2*d0)
            differences[1] = (round(new_voisin_b[1],arrondi) - round(new_voisin_h[1],arrondi))/(2*d0)
        
    elif d0 < d1:
        if exact == 0:
            differences[0] = (round(M_utile[g].x,arrondi) - round(M_utile[g0].x,arrondi))/d0
            differences[1] = (round(M_utile[g].y,arrondi) - round(M_utile[g0].y,arrondi))/d0
        else:
            differences[0] = (round(new_marqueur[0],arrondi) - round(new_voisin_h[0],arrondi))/d0
            differences[1] = (round(new_marqueur[1],arrondi) - round(new_voisin_h[1],arrondi))/d0
    else:
        if exact == 0:
            differences[0] = (round(M_utile[g1].x,arrondi) - round(M_utile[g].x,arrondi))/d1
            differences[1] = (round(M_utile[g1].y,arrondi) - round(M_utile[g].y,arrondi))/d1
        else:
            differences[0] = (round(new_voisin_b[0],arrondi) - round(new_marqueur[0],arrondi))/d1
            differences[1] = (round(new_voisin_b[1],arrondi) - round(new_marqueur[1],arrondi))/d1
    #horizontal
    if d2 == d3:
        if exact == 0:
            differences[2] = (round(M_utile[g3].x,arrondi) - round(M_utile[g2].x,arrondi))/(2*d2)
            differences[3] = (round(M_utile[g3].y,arrondi) - round(M_utile[g2].y,arrondi))/(2*d2)
        else:
            differences[2] = (round(new_voisin_d[0],arrondi) - round(new_voisin_g[0],arrondi))/(2*d2)
            differences[3] = (round(new_voisin_d[1],arrondi) - round(new_voisin_g[1],arrondi))/(2*d2)
    elif d2 < d3:
        if exact == 0:
            differences[2] = round((M_utile[g].x - M_utile[g2].x)/d2,arrondi)
            differences[3] = (round(M_utile[g].y,arrondi) - round(M_utile[g2].y,arrondi))/d2
        else:
            differences[2] = (new_marqueur[0] - new_voisin_g[0])/d2
            differences[3] = (round(new_marqueur[1],arrondi) - round(new_voisin_g[1],arrondi))/d2
    else:
        if exact == 0:
            differences[2] = (round(M_utile[g3].x,arrondi) - round(M_utile[g].x,arrondi))/d3
            differences[3] = (round(M_utile[g3].y,arrondi) - round(M_utile[g].y,arrondi))/d3
        else:
            differences[2] = (round(new_voisin_d[0],arrondi) - round(new_marqueur[0],arrondi))/d3
            differences[3] = (round(new_voisin_d[1],arrondi) - round(new_marqueur[1],arrondi))/d3

    return differences

def flot_forward(marqueur,delta_t,compteur): # for M_ini_utile markers
    new_coord = [marqueur.x,marqueur.y]
    if marqueur.norme() < 1:
        w = 1 - np.sqrt(marqueur.norme())
        t = compteur*delta_t
        x = marqueur.x*np.cos(w*t)-marqueur.y*np.sin(w*t)
        y = marqueur.x*np.sin(w*t)+marqueur.y*np.cos(w*t)
        new_coord = [x,y]
    return new_coord

###############

def distance_euclidienne(x,y,u,v):
    norme = np.sqrt((x-u)*(x-u)+(y-v)*(y-v))
    return norme

def distance_cellule(x,y,u,v,h_ini,quart_plan):
    dist = -1
    #right up :
    if quart_plan == 0:
        if u <= x + h_ini/2:
            dist = np.abs(v - (y + h_ini/2))
        elif v <= y + h_ini/2:
            dist = np.abs(u - (x + h_ini/2))
        else:
            dist = distance_euclidienne(x+h_ini/2,y+h_ini/2,u,v)
    #right down :
    elif quart_plan == 1:
        if u <= x + h_ini / 2:
            dist = np.abs(v - (y - h_ini / 2))
        elif v >= y - h_ini / 2:
            dist = np.abs(u - (x + h_ini / 2))
        else:
            dist = distance_euclidienne(x + h_ini / 2, y - h_ini / 2, u, v)
    #left down :
    elif quart_plan == 2:
        if u >= x - h_ini / 2:
            dist = np.abs(v - (y - h_ini / 2))
        elif v >= y - h_ini / 2:
            dist = np.abs(u - (x - h_ini / 2))
        else:
            dist = distance_euclidienne(x - h_ini / 2, y - h_ini / 2, u, v)
    #left up :
    elif quart_plan == 3:
        if u >= x - h_ini / 2:
            dist = np.abs(v - (y + h_ini / 2))
        elif v <= y + h_ini / 2:
            dist = np.abs(u - (x - h_ini / 2))
        else:
            dist = distance_euclidienne(x - h_ini / 2, y + h_ini / 2, u, v)
    if dist < 0:
        print("error while computing dist")
        exit()
    return dist

def classe_coeff(classement,coeff):
    #descending order
    coeff_copy = []
    for i in range(len(coeff)):
        coeff_copy.append(coeff[i])
    for e in range(3):
        for f in range(4):
            if e < f and coeff_copy[e] < coeff_copy[f]:
                temp1 = classement[e]
                classement[e] = classement[f]
                classement[f] = temp1
                temp2 = coeff_copy[e]
                coeff_copy[e] = coeff_copy[f]
                coeff_copy[f] = temp2

def calcul_marqueur(M_ini_hd,M_ini_bd,M_ini_bg,M_ini_hg,dist_hd,dist_bd,dist_bg,dist_hg,debug):
    a = 0.5
    b = (1-a)/4
    c = 0*(a+b)
    total = 1/dist_hd + 1/dist_bd + 1/dist_bg + 1/dist_hg
    p_hd = b + a/(dist_hd*total)
    p_bd = b + a/(dist_bd*total)
    p_bg = b + a/(dist_bg*total)
    p_hg = b + a/(dist_hg*total)

    max_coeff = max(p_hd,max(p_bd,max(p_bg,p_hg)))
    if max_coeff >= c:

        classement = [0,1,2,3]
        coeff = [p_hd,p_bd,p_bg,p_hg]
        dist = [dist_hd,dist_bd,dist_bg,dist_hg]

        classe_coeff(classement, coeff)

        total2 = 1/dist[classement[1]] + 1/dist[classement[2]] + 1/dist[classement[3]]
        total3 = 1/dist[classement[2]] + 1/dist[classement[3]]
        coeff[classement[1]] = a/(dist[classement[1]]*total) + 3*b/(dist[classement[1]]*total2)

        coeff[classement[2]] = a / (dist[classement[2]] * total) + 3 * b * (
        1 - 1 / (dist[classement[1]] * total2)) * 1 / (dist[classement[2]] * total3)

        coeff[classement[3]] = a / (dist[classement[3]] * total) + 3 * b * (
        1 - 1 / (dist[classement[1]] * total2)) * 1 / (dist[classement[3]] * total3)

        p_hd, p_bd, p_bg, p_hg = coeff

    if debug == -1:
        print("coeff retouchés :")
        print(p_hd)
        print(p_bd)
        print(p_bg)
        print(p_hg)
        print("############")

    x = p_hd * M_ini_hd.x + p_bd * M_ini_bd.x + p_bg * M_ini_bg.x + p_hg * M_ini_hg.x
    y = p_hd * M_ini_hd.y + p_bd * M_ini_bd.y + p_bg * M_ini_bg.y + p_hg * M_ini_hg.y
    marqueur_test = marqueur(x, y, 0, 0, -1, 1)
    return marqueur_test


### Add markers for a given configuration ###
def remaillage(C, h, h_ini, delta_t, compteur, M_ini, M_ini_utile, M_utile, a, b, generation):
    ### Marker at the center of cellule C at time compteur*delta_t ###
    marqueur_test = marqueur(C.x, C.y, C.i, C.j, -1, 1)
    ### We look for the initial position ###
    for e in range(compteur):
        marqueur_test.deplacement_arriere(delta_t)

    h_loc = h_ini / 2
    critere = False
    nb_sub = 1
    while critere == False:
        ### We get the indexes of the 4 closest markers on the originla grid ###
        if marqueur_test.x > 0:
            j = int(marqueur_test.x / h_loc) + int(b / h_loc)
        else:
            j = int(marqueur_test.x / h_loc) + int(b / h_loc) - 1
        if marqueur_test.y > 0:
            i = int(a / h_loc) - int(marqueur_test.y / h_loc)
        else:
            i = int(a / h_loc) - int(marqueur_test.y / h_loc) + 1

        M1 = M_ini[nb_sub][i, j]
        M2 = M_ini[nb_sub][i - 1, j]
        M3 = M_ini[nb_sub][i, j + 1]
        M4 = M_ini[nb_sub][i - 1, j + 1]
        l = marqueur_test.voisins(M1, M2, M3, M4)
        i1 = l[0].i
        j1 = l[0].j

        marqueur_new = marqueur(0, 0, 0, 0, -1, 0)
        marqueur_new.copie(M_ini[nb_sub][i1, j1])
        for e in range(compteur):
            marqueur_new.deplacement(delta_t)
        if C.distance(marqueur_new) <= h_ini / 2:
            critere = True
        else:
            h_loc = h_loc / 2
            nb_sub += 1
            if h_loc < h:
                n = 2 * int(a / h) + 1
                m = 2 * int(b / h) + 1
                M_ini2 = initialisation(a, b, h / 2, 0, -1, None, 0)
                for j in range(m):
                    for i in range(n):
                        if M_ini[-1][i, j].generation >= 0:
                            M_ini2[2 * i, 2 * j].generation = M_ini[-1][i, j].generation
                M_ini.append(M_ini2)
                h = h / 2
                for e in range(len(M_utile)):
                    M_ini_utile[e].i = 2 * M_ini_utile[e].i
                    M_ini_utile[e].j = 2 * M_ini_utile[e].j
                    M_utile[e].i = 2 * M_utile[e].i
                    M_utile[e].j = 2 * M_utile[e].j
    if nb_sub == (len(M_ini) - 1):
        v = j1
        u = i1
    else:
        v = 2 * (len(M_ini) - 1 - nb_sub) * j1
        u = 2 * (len(M_ini) - 1 - nb_sub) * i1
    M_ini[-1][u, v].generation = generation
    marqueur_new_ini = marqueur(0, 0, 0, 0, -1, 0)
    marqueur_new_ini.copie(M_ini[-1][u, v])
    M_ini_utile.append(marqueur_new_ini)
    marqueur_new.generation = generation
    M_utile.append(marqueur_new)
    generation += 1
    variables = []
    variables.append(M_ini)
    variables.append(h)
    variables.append(generation)
    variables.append(marqueur_test)
    return variables

### Add markers for real time ###
def remaillage2(C, h, h_ini, delta_t, compteur, M_ini, M_ini_utile, M_utile, a, b, generation,M_ini_hd,M_ini_bd,M_ini_bg,M_ini_hg,dist_hd,dist_bd,dist_bg,dist_hg):
    h_loc = h_ini / 2
    critere = False
    nb_sub = 1
    debug = 0
    marqueur_test = calcul_marqueur(M_ini_hd, M_ini_bd, M_ini_bg, M_ini_hg, dist_hd, dist_bd, dist_bg, dist_hg,
                                    compteur)

    abs_min = min(M_ini_hd.x, min(M_ini_bd.x, min(M_ini_bg.x, M_ini_hg.x)))
    abs_max = max(M_ini_hd.x, max(M_ini_bd.x, max(M_ini_bg.x, M_ini_hg.x)))
    ord_min = min(M_ini_hd.y, min(M_ini_bd.y, min(M_ini_bg.y, M_ini_hg.y)))
    ord_max = max(M_ini_hd.y, max(M_ini_bd.y, max(M_ini_bg.y, M_ini_hg.y)))
    while critere == False:
        i_start = int(round(a / h_loc)) - int(round(ord_max / h_loc))
        j_start = int(round(abs_min / h_loc)) + int(round(b / h_loc))

        nb_iter_i = 1 + int(round((ord_max-ord_min)/h_loc))
        nb_iter_j = 1 + int(round((abs_max-abs_min)/h_loc))

        for i1 in range(nb_iter_i):
            for j1 in range(nb_iter_j):
                if (i1%2 or j1%2) or nb_sub == 1:
                    marqueur_new = marqueur(0, 0, 0, 0, -1, 0)
                    marqueur_new.copie(M_ini[nb_sub][i_start+i1, j_start+j1])
                    for e in range(compteur):
                        marqueur_new.deplacement(delta_t)
                    if C.distance(marqueur_new) <= (h_ini/2) + 1.E-13:
                        critere = True
                        break
            if critere == True:
                break
        if critere == False:
            h_loc = h_loc / 2
            nb_sub += 1
            if h_loc < h:
                n = 2 * int(a / h) + 1
                m = 2 * int(b / h) + 1
                M_ini2 = initialisation(a, b, h / 2, 0, -1, None, 0)
                for j in range(m):
                    for i in range(n):
                        if M_ini[-1][i, j].generation >= 0:
                            M_ini2[2 * i, 2 * j].generation = M_ini[-1][i, j].generation
                M_ini.append(M_ini2)
                h = h / 2
                for e in range(len(M_utile)):
                    M_ini_utile[e].i = 2 * M_ini_utile[e].i
                    M_ini_utile[e].j = 2 * M_ini_utile[e].j
                    M_utile[e].i = 2 * M_utile[e].i
                    M_utile[e].j = 2 * M_utile[e].j
    if nb_sub == (len(M_ini) - 1):
        v = j_start+j1
        u = i_start+i1
    else:
        v = 2 * (len(M_ini) - 1 - nb_sub) * (j_start+j1)
        u = 2 * (len(M_ini) - 1 - nb_sub) * (i_start+i1)
    M_ini[-1][u, v].generation = generation
    marqueur_new_ini = marqueur(0, 0, 0, 0, -1, 0)
    marqueur_new_ini.copie(M_ini[-1][u, v])
    M_ini_utile.append(marqueur_new_ini)
    marqueur_new.generation = generation
    M_utile.append(marqueur_new)
    generation += 1
    variables = []
    variables.append(M_ini)
    variables.append(h)
    variables.append(generation)
    variables.append(marqueur_test)
    return variables