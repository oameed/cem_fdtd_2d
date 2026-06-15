###################################################
### FINITE DIFFERENCE TIME DOMAIN (FDTD) 1D TEz ###
### SOLVER WITH:                                ###
###     POINT SOURCE EXCITATION                 ###
###     SPLIT-PML ABC                           ###
### by: OAMEED NOAKOASTEEN                      ###
###################################################

import os
import numpy as np

def rJSON(filename):
    import json
    with open(filename, 'r') as file:
        x = json.load(file)
    return x

def get_wavelength_min(params):
    def get_wavelength_min_gaussian(params,c0):
        fmax = params[1]['freq']
        return c0/fmax
    def get_wavelength_min_cmg(params,c0):
        bw   = params[1]['bw (pct)']/100
        fmax = params[1]['freq']*(1 + 0.5*bw)
        return c0/fmax
    options             = {'gaussian': get_wavelength_min_gaussian,
                           'cmg'     : get_wavelength_min_cmg      }
    wavelength_min_func = options[params[1]['pulse']]
    c0                  = 1/np.sqrt(params[2][1]*params[2][0])
    return wavelength_min_func(params,c0)

def get_grid_spec(params):
    # grid_spec[0]: N
    # grid_spec[1]: delta
    # grid_spec[2]: dt
    c0                      = 1/np.sqrt(params[2][1]*params[2][0])
    lx                      = 1
    wavelength_min          = get_wavelength_min(params)
    delta                   = wavelength_min/params[0][0]
    N                       = int(np.floor(lx/delta))
    dt                      = params[2][2]*(delta/c0)
    return [N,delta,dt] 

def get_pulse(grid_spec,params):
    def get_pulse_gaussian(time,params):
        taw  = np.sqrt(2.3)/(np.pi*params[1]['freq'])
        t0   = np.sqrt(20)*taw
        func = lambda t: np.exp(-np.power((t-t0)/taw,2))
        return [func(t) for t in time]
    def get_pulse_cmg(time,params):
        bw   = params[1]['bw (pct)']/100
        df   = bw*params[1]['freq']
        taw  = np.sqrt(2.3)/(np.pi*(df/2))
        t0   = np.sqrt(20)*taw
        func = lambda t: np.cos(2*np.pi*params[1]['freq']*(t-t0))*np.exp(-np.power((t-t0)/taw,2))
        return [func(t) for t in time]
    options_type    = {'gaussian'   : get_pulse_gaussian,
                       'cmg'        : get_pulse_cmg      }
    options_loc     = {'cavity'     : 0.5 ,
                       'propagation': 0.32 }
    pulse_type_func = options_type[params[1]['pulse'       ]]
    loc_factor      = options_loc [params[1]['problem_type']]
    pulse_loc       = int(loc_factor*grid_spec[0])
    time            = [grid_spec[2]*n for n in range(params[0][1])]
    pulse           = pulse_type_func(time,params)
    return {'val':pulse, 'loc':pulse_loc}

def get_pml_region(MPey,MPhz,grid_spec,params):
    def get_pml_region_cavity(MPey,MPhz,grid_spec,params):
        pass
    def get_pml_region_propagation(MPey,MPhz,grid_spec,params):
        size_factor           = 0.25
        R0                    = 1e-9
        P                     = 2
        thickness             = int(size_factor*grid_spec[0])
        sigmamax              = -1*np.sqrt(params[2][0]/params[2][1])*(((P+1)*np.log(R0))/(2*thickness*grid_spec[1]))
        start = 0
        end   = thickness
        for x in [x for x in range(start,end)]:
            MPey[x,1] = sigmamax*np.power(((end-1) - x)/((end-1) - start),P)
        start = grid_spec[0] - thickness
        end   = grid_spec[0]
        for x in [x for x in range(start,end)]:
            MPey[x,1] = sigmamax*np.power((x - start)/((end-1) - start),P)
        MPhz[:,1] = (params[2][1]/params[2][0])*MPey[:,1]
    
    options  = {'cavity'     : get_pml_region_cavity     ,
                'propagation': get_pml_region_propagation }
    pml_func = options[params[1]['problem_type']]
    pml_func(MPey,MPhz,grid_spec,params)

def get_field_coefficients(Cey,Chz,grid_spec,params):
    # MPey[:,0]: PERMATIVITTY
    # MPey[:,1]: ELECTRIC CONDUCTIVITY
    # MPhz[:,0]: PERMEABILITY
    # MPhz[:,1]: MAGNETIC CONDUCTIVITY
    # Cey [:,0]: Ceye
    # Cey [:,1]: Ceyh
    # Cey [:,2]: Ceyj
    # Chz [:,0]: Chzh
    # Chz [:,2]: Chzey
    # Chz [:,3]: Chzm
    MPey      =  np.zeros((grid_spec[0],2))
    MPhz      =  np.zeros((grid_spec[0],2))
    MPey[:,0] =  params[2][0]
    MPhz[:,0] =  params[2][1]
    get_pml_region(MPey,MPhz,grid_spec,params)
    Cey [:,0] =  (2*MPey[:,0] - grid_spec[2]*MPey[:,1])/( 2*MPey[:,0] + grid_spec[2]*MPey[:,1]              ) 
    Cey [:,1] = -(2*            grid_spec[2]          )/((2*MPey[:,0] + grid_spec[2]*MPey[:,1])*grid_spec[1])
    Cey [:,2] = -(2*            grid_spec[2]          )/( 2*MPey[:,0] + grid_spec[2]*MPey[:,1]              )
    Chz [:,0] =  (2*MPhz[:,0] - grid_spec[2]*MPhz[:,1])/( 2*MPhz[:,0] + grid_spec[2]*MPhz[:,1]              )
    Chz [:,1] = -(2*            grid_spec[2]          )/((2*MPhz[:,0] + grid_spec[2]*MPhz[:,1])*grid_spec[1])
    Chz [:,2] = -(2*            grid_spec[2]          )/((2*MPhz[:,0] + grid_spec[2]*MPhz[:,1])             )

def extract_view_and_calculate_output(ey,hz,grid_spec,params):
    def extract_view_cavity(ey,hz,grid_spec):
        return ey.copy(), hz.copy()
    def extract_view_propagation(ey,hz,grid_spec):
        size_factor     = 0.28
        thickness       = int(size_factor*grid_spec[0])
        index_start     = thickness + 1
        index_end       = grid_spec[0] - thickness
        ey_copy         = ey.copy()
        ey_copy         = ey_copy[index_start:index_end]
        hz_copy         = hz.copy()
        hz_copy         = hz_copy[index_start:index_end]
        return ey_copy, hz_copy
    def calculate_output(ey,hz):
        ey_copy = ey
        return ey_copy
    options           = {'cavity'     : extract_view_cavity     ,
                         'propagation': extract_view_propagation }
    extract_view_func = options[params[1]['problem_type']]
    ey_copy, hz_copy  = extract_view_func(ey,hz,grid_spec)
    return calculate_output(ey_copy,hz_copy)

def visualizer(exportarray,params,paths):
    def update(frame_number,exportarray,imgobj,txtobj):
        imgobj.set_ydata(exportarray[frame_number])
        txtobj.set_text ('{}'.format(frame_number))
        return [imgobj,txtobj]
    import matplotlib.pyplot    as plt
    import matplotlib.animation as animation
    size_base  = 4
    size_ratio = 0.25
    fps        = 10
    dpi        = 90
    filename   = os.path.join(paths[0],'data' + '.gif')
    arraymax   = max([np.amax(np.abs(x)) for x in exportarray]) 
    fig, ax    = plt.subplots(figsize=(size_base,size_ratio*size_base))
    ax.set_xlim(0,len(exportarray[0])-1)
    ax.set_ylim(-arraymax,arraymax)
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, wspace = 0, hspace = 0)
    ax.axis('off')
    imgobj     = ax.plot(exportarray[0])[0]
    txtobj     = ax.text(0.05,0.05,'', transform = ax.transAxes, color = 'k', fontsize = 8, fontweight = 'bold')
    updatefunc = lambda n: update(n,exportarray,imgobj,txtobj)
    animobj    = animation.FuncAnimation(fig = fig, func = updatefunc, frames = len(exportarray))
    animobj.save(filename, writer = 'pillow', dpi = dpi, fps = fps)
    plt.close()

def initialize_run():
    # params[0][0]: number of cells per wavelength
    # params[0][1]: number of time steps
    # params[0][2]: problem type
    # params[1][*]: dictionay containing structure information
    # params[2][0]: FREE SPACE PERMITTIVITY
    # params[2][1]: FREE SPACE PERMEABILITY
    # params[2][2]: CFL FACTOR
    import argparse
    import shutil
    parser     = argparse.ArgumentParser()
    parser.add_argument('-p', type = str, required = True)
    parser.add_argument('-v', type = str, required = True)
    parser.add_argument('-n', type = int, default  = 40  )
    parser.add_argument('-t', type = int, default  = 1000)
    args       = parser.parse_args()
    structure  = rJSON(args.p + '.json')
    params     = [[args.n,args.t,args.p]       ,
                  structure                    ,
                  [8.85e-12, 4*np.pi*1e-7, 0.5] ]
    paths      = [os.path.join('..','experiments',args.v)]
    shutil.rmtree(paths[0], ignore_errors = True)
    for path in paths:
        os.makedirs(path)
    return params, paths

def main():
    params, paths = initialize_run()
    exportarray   = []
    
    grid_spec     = get_grid_spec(params)
    
    Cey           = np.zeros((grid_spec[0],3))
    Chz           = np.zeros((grid_spec[0],3))
    get_field_coefficients(Cey,Chz,grid_spec,params)
    
    Ey            = np.zeros((grid_spec[0],))
    Hz            = np.zeros((grid_spec[0],))
    
    pulse         = get_pulse(grid_spec,params)
    
    print("Begin time updates ...")
    for n in range(params[0][1]):
        
        i     = pulse['loc']
        Hz[i] = Hz[i] + Chz[i,2]*pulse['val'][n]
        
        for i in range(  grid_spec[0]-1):
            Hz[i] = Chz[i,0]*Hz[i] + Chz[i,1]*(Ey[i+1] - Ey[i  ])

        for i in range(1,grid_spec[0]-1):
            Ey[i] = Cey[i,0]*Ey[i] + Cey[i,1]*(Hz[i  ] - Hz[i-1])

        exportarray.append(extract_view_and_calculate_output(Ey,Hz,grid_spec,params))
    visualizer(exportarray,params,paths)
    print("Finished!")

if __name__ == "__main__":
    main()
