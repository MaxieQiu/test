
# [file name]: Cubic_Spline.py
# [file content begin]
import numpy as np
import matplotlib.pyplot as plt

def Cubic_Spline_interpolation(coordinates):
    x, y = zip(*sorted(coordinates, key=lambda p: p[0]))
    x, y = np.array(x), np.array(y)
    n = len(x) - 1
    h = np.diff(x)
    delta = np.diff(y) / h
    
    
    main_diag = 2 * (h[:-1] + h[1:])
    upper_diag = h[1:-1] if n > 2 else h[1:]
    lower_diag = h[1:-1] if n > 2 else h[:-1]
    
    
    c = np.zeros(n+1)
    for i in range(1, n-1):
        lower = lower_diag[i-1]/main_diag[i-1]
        main_diag[i] -= lower * upper_diag[i-1]
        delta[i] -= lower * delta[i-1]
    c[-2] = delta[-1]/main_diag[-1]
    for i in range(n-3, -1, -1):
        c[i+1] = (delta[i] - upper_diag[i] * c[i+2])/main_diag[i]
    
   
    b = delta - h*(2*c[:-1] + c[1:])/6
    d = (c[1:] - c[:-1])/(6*h)
    
    return lambda z: [y[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])**2/2 + d[i]*(t-x[i])**3/6 
                     for t in z 
                     for i in [np.clip(np.searchsorted(x, t, side='right')-1, 0, n-1)] 
                     if x[i] <= t <= x[i+1]]

def tests_set_of_points(coordinates):
    x_i = np.array([p[0] for p in coordinates])
    y_i = np.array([p[1] for p in coordinates])
    g = Cubic_Spline_interpolation(coordinates)
    x = np.linspace(min(x_i), max(x_i), 100)
    y = []
    for t in x:
        try:
            y.append(g([t])[0])
        except:
            y.append(np.nan)
    plt.plot(x, y, color='b', label='approximation')
    plt.scatter(x_i, y_i, s=100, color='r', label='points')
    plt.legend()
    plt.show()

coordinates = [(-5,-3),(-3,-1),(-1,0),(0,2),(1,1),(2,4),(3,9),(4,8),(5,5)]
tests_set_of_points(coordinates)

