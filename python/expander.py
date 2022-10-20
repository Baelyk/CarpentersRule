import numpy as np
from statistics import mean
import matplotlib.pyplot as plt
import math, scipy
import cvxpy as cp

def max_chord_length(coords, return_list=False):
    lengths = [np.linalg.norm(coords[i]-coords[i+1]) for i in range(len(coords)-1)]
    if return_list:
        return lengths, sum(lengths)
    return sum(lengths)

def compute_agg_lengths(lengths):
    # a list of aggregated edge length
    return [sum(lengths[:i+1]) for i in range(len(lengths))]

def get_edge_ind(agg_lengths, a):
    if a>=agg_lengths[-1]:
        return None
    for n in range(len(agg_lengths)):
        if agg_lengths[n]>a:
            return n

def points_on_adjacent_edge(agg_lengths, a, b):
    try:
        if abs(get_edge_ind(agg_lengths, a) - get_edge_ind(agg_lengths, b)) == 1:
            return True
        else:
            return False
    except:
        return True


def gamma(coords, b, agg_lengths=None):
    # get or calculate a list of edge length and total edge length
    if agg_lengths==None:
        lengths, max_l = max_chord_length(coords, return_list=True)
        agg_lengths = [sum(lengths[:i+1]) for i in range(len(lengths))]
    # return infinity if input longer than total edge length
    if b>=agg_lengths[-1]: # 0<=a<b<l
        return float('inf') #"b too big!!"
    n = get_edge_ind(agg_lengths, b)
    if n==None:
        return None
#     print(n)
    if n==0:
        portion = 1-(b/agg_lengths[0])
    else:
        portion = 1-((b-agg_lengths[n-1])/(agg_lengths[n]-agg_lengths[n-1]))
    x=coords[n]
    y=coords[n+1]
    return portion*x+(1-portion)*y

def polygon(n):
    points = np.random.rand(n,2)
    for i in range(n-1):
        for j in range(i+1,n):
            if points[i].all()==points[j].all():
                points[j] = np.random.rand(1,2)
    # reordering vertices
#     print(points)
    x_c = points[:,0]
    y_c = points[:,1]
    mean_x = mean(x_c)
    mean_y = mean(y_c)
    angles = []
    for i in range(len(points)):
        angle = math.atan2( (y_c[i]-mean_y),(x_c[i]-mean_x) )
        angles.append(angle)
    sort_index = np.argsort(np.array(angles))
    new_points = points[sort_index]
    return new_points
# coords = polygon(200)

def spiral(e, display=False, linewidth=1):
    # length of coords = 2/e+2
    coords = [[0,0], [1,0], [1,1], [0,1]]
    ns = [None]*int(2/e-2)
    coords.extend(ns)
    for i in range(4,len(coords)):
#         if i>=2/e+2:
#             coords=coords[:i]
#             plt.plot([i[0] for i in coords],[i[1] for i in coords], '-o')
#             break
        if i%4==0:
            coords[i] = [coords[i-1][0],coords[i-4][1]+e]
        elif i%4==1:
            coords[i] = [coords[i-4][0]-e,coords[i-1][1]]
        elif i%4==2:
            coords[i] = [coords[i-1][0],coords[i-4][1]-e]
        else:
            coords[i] = [coords[i-4][0]+e,coords[i-1][1]]
    if display==True:
        fig, ax = plt.subplots()
        ax.axis('off')
        fig.patch.set_visible(False)
        ax.plot([i[0] for i in coords],[i[1] for i in coords], '-o', linewidth=linewidth)
    return np.array([np.array(i) for i in coords])

def generate_constraints(coords, lb_=0, ub_=float('inf')):
    # A is a m-by-2n matrix, where m is the number of constraints and n=|V| 
    A=[]
    n = len(coords)
    m = int(n*(n-1)/2)
    lb=[0]*m
    ub=[0]*m
    ind_m=0
    for i in range(n-1):
        for j in range(i+1,n):
            # indices of vix, viy, vjx, vjy
            ind_vix = i*2
            ind_viy = ind_vix+1
            ind_vjx = j*2
            ind_vjy = ind_vjx+1
            # pix, piy, pjx, pjy
            pix = coords[i][0]
            piy = coords[i][1]
            pjx = coords[j][0]
            pjy = coords[j][1]
            # row
            row=[0]*(2*n)
            row[ind_vix] = pix-pjx
            row[ind_viy] = piy-pjy
            row[ind_vjx] = pjx-pix
            row[ind_vjy] = pjy-piy
            if j-i!=1:
                lb[ind_m] = lb_ #0#e-10 #np.linalg.norm(coords[j]-coords[i])+2**(-10)
                ub[ind_m] = ub_ #float('inf')
            ind_m+=1
            A.append(row)
    return np.array(A), lb, ub

def obj_fun_generator(coords,first=False,second=False,both=False):
    coords=np.array([np.array(i) for i in coords])
    x = cp.Variable(2*len(coords))
    p=coords
    part1 = cp.sum_squares(x)
    # part1 = sum([i**2 for i in x])
    # part1=0
    part2 = 0
    for i in range(len(p) - 2):
        for j in range(i+2, len(p)):
    #             print(i, j)
            vi = np.array([x[i*2], x[i*2+1]])
            vj = np.array([x[j*2], x[j*2+1]])
            part2a = (vi-vj)[0]*(p[i]-p[j])[0] + (vi-vj)[1]*(p[i]-p[j])[1]
            den = part2a - cp.norm(p[j]-p[i]) # math.sqrt(sum([k**2 for k in p[j]-p[i]]))
            part2+=den**(-1)
    if first:
        return part1
    elif second:
        return part2
    elif both:
        return part1+part2

def expand(coords, x, obj_fun, constraints, iteration=40, interval=10, step_size=1e-3, method='CPLEX'):
    coords=np.array([np.array(i) for i in coords])
    plt.plot([i[0] for i in coords],[i[1] for i in coords], '-o')
    q=0
#     x = cp.Variable(2*len(coords))
    for _ in range(40):
        problem = cp.Problem(cp.Minimize(obj_fun), constraints)
        problem.solve(solver= method, verbose=False)
        for i in range(len(x.value)):
            coords[int(i/2)][i%2]=coords[int(i/2)][i%2]+(x.value[i]*step_size)
        # plot
        q+=1
        if q%interval==0:
            plt.plot([i[0] for i in coords],[i[1] for i in coords], '-o')
    # solver= "CPLEX", 