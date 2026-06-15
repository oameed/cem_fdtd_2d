###################################################
### FINITE DIFFERENCE TIME DOMAIN (FDTD) 2D TEz ###
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
    # grid_spec[0]: Nx
    # grid_spec[1]: Ny
    # grid_spec[2]: delta
    # grid_spec[3]: dt
    options_ly_scale_factor = {'cavity'     : params[1]['ratio'],
                               'propagation': 1                  }
    ly_scale_factor         = options_ly_scale_factor[params[1]['problem_type']]
    c0                      = 1/np.sqrt(params[2][1]*params[2][0])
    lx                      = 1
    ly                      = ly_scale_factor*lx
    wavelength_min          = get_wavelength_min(params)
    delta                   = wavelength_min/params[0][0]
    Nx                      = int(np.floor(lx/delta))
    Ny                      = int(np.floor(ly/delta))
    dt                      = params[2][2]*(delta/c0)
    return [Nx,Ny,delta,dt] 

def get_update_indeces(grid_spec):
    def get_indeces_ex(grid_spec):
        indeces = []
        for     i in range(1,grid_spec[0]-1):
            for j in range(  grid_spec[1]-1):
                indeces.append((i,j))
        return indeces
    def get_indeces_ey(grid_spec):
        indeces = []
        for     i in range(  grid_spec[0]-1):
            for j in range(1,grid_spec[1]-1):
                indeces.append((i,j))
        return indeces
    def get_indeces_hz(grid_spec):
        size_factor     = 0.25
        thickness_x     = int(size_factor*grid_spec[0])
        thickness_y     = int(size_factor*grid_spec[1])
        indeces         = []
        indeces_inside  = []
        for     i in range(grid_spec[0]-1):
            for j in range(grid_spec[1]-1):
                indeces.append((i,j))
        for     i in range(thickness_x + 1, grid_spec[0] - thickness_x):
            for j in range(thickness_y + 1, grid_spec[1] - thickness_y):
                indeces_inside.append((i,j))
        indeces_outside = [x for x in indeces if not x in indeces_inside]
        return indeces_inside, indeces_outside 
    indeces_ex                            = get_indeces_ex(grid_spec)
    indeces_ey                            = get_indeces_ey(grid_spec)
    indeces_hz_inside, indeces_hz_outside = get_indeces_hz(grid_spec)
    return {'ex':indeces_ex, 'ey':indeces_ey, 'hz_inside':indeces_hz_inside, 'hz_outside':indeces_hz_outside}

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
    def get_pulse_loc(loc_factor,grid_spec):
        locx       = int(loc_factor*grid_spec[0])
        locy       = int(loc_factor*grid_spec[1])
        #locy       = grid_spec[1] - locy
        return [locx,locy]
    options_type    = {'gaussian'   : get_pulse_gaussian,
                       'cmg'        : get_pulse_cmg      }
    options_loc     = {'cavity'     : 0.5 ,
                       'propagation': 0.32 }
    pulse_type_func = options_type[params[1]['pulse'       ]]
    loc_factor      = options_loc [params[1]['problem_type']]
    time            = [grid_spec[3]*n for n in range(params[0][1])]
    pulse           = pulse_type_func(time,params)
    pulse_loc       = get_pulse_loc(loc_factor,grid_spec)
    return {'val':pulse, 'loc':pulse_loc}

def get_objects(MPex,MPey,MPhz,grid_spec,params):
    def get_objects_cavity(MPex,MPey,MPhz,grid_spec,params): 
        def get_objects_noobjects(MPex,MPey,MPhz,grid_spec,params): 
            pass
        def get_objects_circle(MPex,MPey,MPhz,grid_spec,params):
            center_i = int(grid_spec[0]/2)
            center_j = int(grid_spec[1]/2)
            radius   = params[1]['radius (m)']/grid_spec[2]
            for     i in range(grid_spec[0]):
                for j in range(grid_spec[1]):
                    if np.sqrt(np.power(i-center_i,2)+np.power(j-center_j,2)) > radius:
                        MPex[i,j,1] = params[1]['sigma']
                        MPey[i,j,1] = params[1]['sigma']
        options            = {'rectangle': get_objects_noobjects,
                              'circle'   : get_objects_circle    }
        object_cavity_func = options[params[1]['name']]
        object_cavity_func(MPex,MPey,MPhz,grid_spec,params)
    def get_objects_propagation(MPex,MPey,MPhz,grid_spec,params):
        def get_objects_noobjects(MPex,MPey,MPhz,grid_spec,params):
            pass
        def get_objects_rectanglePEC(MPex,MPey,MPhz,grid_spec,params):
            center_i = int(grid_spec[0]/2)
            center_j = int(grid_spec[1]/2)
            side     = params[1]['side (m)']/grid_spec[2]
            for     i in range(grid_spec[0]):
                for j in range(grid_spec[1]):
                    if center_i - side/2 < i < center_i + side/2 and center_j - side/2 < j < center_j + side/2:
                        MPex[i,j,1] = params[1]['sigma']
                        MPey[i,j,1] = params[1]['sigma']
        def get_objects_rectangleDE(MPex,MPey,MPhz,grid_spec,params):
            center_i = int(grid_spec[0]/2)
            center_j = int(grid_spec[1]/2)
            side     = params[1]['side (m)']/grid_spec[2]
            for     i in range(grid_spec[0]):
                for j in range(grid_spec[1]):
                    if center_i - side/2 < i < center_i + side/2 and center_j - side/2 < j < center_j + side/2:
                        MPex[i,j,0] = MPex[i,j,0]*params[1]['epsr']
                        MPey[i,j,0] = MPey[i,j,0]*params[1]['epsr']
        def get_objects_circlePEC(MPex,MPey,MPhz,grid_spec,params):
            center_i = int(grid_spec[0]/2)
            center_j = int(grid_spec[1]/2)
            radius   = params[1]['radius (m)']/grid_spec[2]
            for     i in range(grid_spec[0]):
                for j in range(grid_spec[1]):
                    if np.sqrt(np.power(i-center_i,2)+np.power(j-center_j,2)) < radius:
                        MPex[i,j,1] = params[1]['sigma']
                        MPey[i,j,1] = params[1]['sigma']
        def get_objects_circleDE(MPex,MPey,MPhz,grid_spec,params):
            center_i = int(grid_spec[0]/2)
            center_j = int(grid_spec[1]/2)
            radius   = params[1]['radius (m)']/grid_spec[2]
            for     i in range(grid_spec[0]):
                for j in range(grid_spec[1]):
                    if np.sqrt(np.power(i-center_i,2)+np.power(j-center_j,2)) < radius:
                        MPex[i,j,0] = MPex[i,j,0]*params[1]['epsr']
                        MPey[i,j,0] = MPey[i,j,0]*params[1]['epsr']
        options                 = {'noobjects'   : get_objects_noobjects   ,
                                   'rectanglePEC': get_objects_rectanglePEC,
                                   'rectangleDE' : get_objects_rectangleDE ,
                                   'circlePEC'   : get_objects_circlePEC   ,
                                   'circleDE'    : get_objects_circleDE     }
        object_propagation_func = options[params[1]['name']]
        object_propagation_func(MPex,MPey,MPhz,grid_spec,params)

    options     = {'cavity'     : get_objects_cavity     ,
                   'propagation': get_objects_propagation }
    object_func = options[params[1]['problem_type']]
    object_func(MPex,MPey,MPhz,grid_spec,params)

def get_pml_region(MPex,MPey,MPhz,grid_spec,params):
    def get_pml_region_cavity(MPex,MPey,MPhz,grid_spec,params):
        pass
    def get_pml_region_propagation(MPex,MPey,MPhz,grid_spec,params):
        size_factor           = 0.25
        R0                    = 1e-9
        P                     = 2
        thickness             = int(size_factor*grid_spec[0])
        sigmamax              = -1*np.sqrt(params[2][0]/params[2][1])*(((P+1)*np.log(R0))/(2*thickness*grid_spec[2]))
        start = 0
        end   = thickness
        for x in [x for x in range(start,end)]:
            MPey[:,x,1] = sigmamax*np.power(((end-1) - x)/((end-1) - start),P)
        start = grid_spec[0] - thickness
        end   = grid_spec[0]
        for x in [x for x in range(start,end)]:
            MPey[:,x,1] = sigmamax*np.power((x - start)/((end-1) - start),P)
        start = 0
        end   = thickness
        for x in [x for x in range(start,end)]:
            MPex[x,:,1] = sigmamax*np.power(((end-1) - x)/((end-1) - start),P)
        start = grid_spec[1] - thickness
        end   = grid_spec[1]
        for x in [x for x in range(start,end)]:
            MPex[x,:,1] = sigmamax*np.power((x - start)/((end-1) - start),P)
        MPhz[:,:,2] = (params[2][1]/params[2][0])*MPey[:,:,1]
        MPhz[:,:,3] = (params[2][1]/params[2][0])*MPex[:,:,1]
    
    options  = {'cavity'     : get_pml_region_cavity     ,
                'propagation': get_pml_region_propagation }
    pml_func = options[params[1]['problem_type']]
    pml_func(MPex,MPey,MPhz,grid_spec,params)

def get_field_coefficients(Cex,Cey,Chz,grid_spec,params):
    # MPex[:,:,0]: PERMATIVITTY
    # MPex[:,:,1]: ELECTRIC CONDUCTIVITY IN NON-PML REGION / IN PML REGION: sigma_e_y
    # MPey[:,:,0]: PERMATIVITTY
    # MPey[:,:,0]: ELECTRIC CONDUCTIVITY IN NON-PML REGION / IN PML REGION: sigma_e_x
    # MPhz[:,:,0]: PERMEABILITY
    # MPhz[:,:,1]: MAGNETIC CONDUCTIVITY IN NON-PML REGION
    # MPhz[:,:,2]: MAGNETIC CONDUCTIVITY IN     PML REGION: sigma_m_x 
    # MPhz[:,:,3]: MAGNETIC CONDUCTIVITY IN     PML REGION: sigma_m_y
    # Cex [:,:,0]: Cexe
    # Cex [:,:,1]: Cexh
    # Cex [:,:,2]: Cexj
    # Cey [:,:,0]: Ceye
    # Cey [:,:,1]: Ceyh
    # Cey [:,:,2]: Ceyj
    # Chz [:,:,0]: Chzh   (NON-PML REGION)
    # Chz [:,:,1]: Chzex  (NON-PML REGION)
    # Chz [:,:,2]: Chzey  (NON-PML REGION)
    # Chz [:,:,3]: Chzm   (NON-PML REGION)
    # Chz [:,:,4]: Chzhzx (    PML REGION)
    # Chz [:,:,5]: Chzey  (    PML REGION)
    # Chz [:,:,6]: Chzhzy (    PML REGION)
    # Chz [:,:,7]: Chzex  (    PML REGION)
    MPex         = np.zeros((grid_spec[0],grid_spec[1],2))
    MPey         = np.zeros((grid_spec[0],grid_spec[1],2))
    MPhz         = np.zeros((grid_spec[0],grid_spec[1],4))
    MPex[:,:,0]  = params[2][0]
    MPey[:,:,0]  = params[2][0]
    MPhz[:,:,0]  = params[2][1]
    get_pml_region(MPex,MPey,MPhz,grid_spec,params)
    get_objects   (MPex,MPey,MPhz,grid_spec,params)
    Cex [:,:,0]  =  (2*MPex[:,:,0] - grid_spec[3]*MPex[:,:,1])/( 2*MPex[:,:,0] + grid_spec[3]*MPex[:,:,1]              )
    Cex [:,:,1]  =  (2*              grid_spec[3]            )/((2*MPex[:,:,0] + grid_spec[3]*MPex[:,:,1])*grid_spec[2])
    Cex [:,:,2]  = -(2*              grid_spec[3]            )/( 2*MPex[:,:,0] + grid_spec[3]*MPex[:,:,1]              )
    Cey [:,:,0]  =  (2*MPey[:,:,0] - grid_spec[3]*MPey[:,:,1])/( 2*MPey[:,:,0] + grid_spec[3]*MPey[:,:,1]              ) 
    Cey [:,:,1]  = -(2*              grid_spec[3]            )/((2*MPey[:,:,0] + grid_spec[3]*MPey[:,:,1])*grid_spec[2])
    Cey [:,:,2]  = -(2*              grid_spec[3]            )/( 2*MPey[:,:,0] + grid_spec[3]*MPey[:,:,1]              )
    Chz [:,:,0]  =  (2*MPhz[:,:,0] - grid_spec[3]*MPhz[:,:,1])/( 2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,1]              )
    Chz [:,:,1]  =  (2*              grid_spec[3]            )/((2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,1])*grid_spec[2])
    Chz [:,:,2]  = -(2*              grid_spec[3]            )/((2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,1])*grid_spec[2])
    Chz [:,:,3]  = -(2*              grid_spec[3]            )/ (2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,1]              )
    Chz [:,:,4]  =  (2*MPhz[:,:,0] - grid_spec[3]*MPhz[:,:,2])/( 2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,2]              )
    Chz [:,:,5]  = -(2*              grid_spec[3]            )/((2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,2])*grid_spec[2])
    Chz [:,:,6]  =  (2*MPhz[:,:,0] - grid_spec[3]*MPhz[:,:,3])/( 2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,3]              )
    Chz [:,:,7]  =  (2*              grid_spec[3]            )/((2*MPhz[:,:,0] + grid_spec[3]*MPhz[:,:,3])*grid_spec[2])

def extract_view_and_calculate_output(ex,ey,hz,grid_spec,params):
    def extract_view_cavity(ex,ey,hz,grid_spec):
        ex_copy = np.transpose(ex,[1,0])
        ey_copy = np.transpose(ey,[1,0])
        hz_copy = np.transpose(hz,[1,0])
        return ex_copy, ey_copy, hz_copy
    def extract_view_propagation(ex,ey,hz,grid_spec):
        size_factor     = 0.28
        thickness_x     = int(size_factor*grid_spec[0])
        thickness_y     = int(size_factor*grid_spec[1])
        index_start_x   = thickness_x + 1
        index_end_x     = grid_spec[0] - thickness_x
        index_start_y   = thickness_y + 1
        index_end_y     = grid_spec[1] - thickness_y
        ex_copy         = ex[index_start_x:index_end_x, index_start_y:index_end_y]
        ex_copy         = np.transpose(ex_copy,[1,0])
        ey_copy         = ey[index_start_x:index_end_x, index_start_y:index_end_y]
        ey_copy         = np.transpose(ey_copy,[1,0])
        hz_copy         = hz[index_start_x:index_end_x, index_start_y:index_end_y]
        hz_copy         = np.transpose(hz_copy,[1,0])
        return ex_copy, ey_copy, hz_copy
    def calculate_output(ex,ey,hz):
        return 0.5*np.abs(hz)*np.sqrt(np.power(np.abs(ex),2)+np.power(np.abs(ey),2))
    options                   = {'cavity'     : extract_view_cavity     ,
                                 'propagation': extract_view_propagation }
    extract_view_func         = options[params[1]['problem_type']]
    ex_copy, ey_copy, hz_copy = extract_view_func(ex,ey,hz,grid_spec)
    return calculate_output(ex_copy,ey_copy,hz_copy)

def visualizer(exportarray,params,paths):
    def get_color_maximum(exportarray,params):
        level_db         = params[0][2] 
        arraymax         = max([np.amax(x) for x in exportarray]) 
        arraymax_db      = 10*np.log10(arraymax)
        color_maximum_db = arraymax_db - level_db
        return np.power(10,color_maximum_db/10)
    def update(frame_number,exportarray,imgobj,txtobj):
        imgobj.set_array(exportarray[frame_number])
        txtobj.set_text ('{}'.format(frame_number))
        return [imgobj,txtobj]
    import matplotlib.pyplot    as plt
    import matplotlib.animation as animation
    options_size_ratio = {'cavity'     : params[1]['ratio'],
                          'propagation': 1                  }
    size_base  = 4
    size_ratio = options_size_ratio[params[1]['problem_type']]
    cmap       = 'jet'
    fps        = 10
    dpi        = 90
    filename   = os.path.join(paths[0],'data' + '.gif')
    vmax       = get_color_maximum(exportarray,params) 
    fig, ax    = plt.subplots(figsize=(size_base,size_ratio*size_base))
    fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, wspace = 0, hspace = 0)
    ax.axis('off')
    imgobj     = ax.imshow(exportarray[0], cmap = cmap, vmin = 0, vmax = vmax, interpolation = 'bilinear')
    txtobj     = ax.text(0.05,0.05,'', transform = ax.transAxes, color = 'w', fontsize = 8, fontweight = 'bold')
    updatefunc = lambda n: update(n,exportarray,imgobj,txtobj)
    animobj    = animation.FuncAnimation(fig = fig, func = updatefunc, frames = len(exportarray))
    animobj.save(filename, writer = 'pillow', dpi = dpi, fps = fps)
    plt.close()

def initialize_run():
    # params[0][0]: number of cells per wavelength
    # params[0][1]: number of time steps
    # params[0][2]: visualization color scale level
    # params[0][3]: problem type
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
    parser.add_argument('-t', type = int, default  = 500 )
    parser.add_argument('-l', type = int, default  = 25  )
    args       = parser.parse_args()
    structure  = rJSON(args.p + '.json')
    params     = [[args.n,args.t,args.l,args.p],
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
    
    indeces       = get_update_indeces(grid_spec)
    
    Cex           = np.zeros((grid_spec[0],grid_spec[1],3))
    Cey           = np.zeros((grid_spec[0],grid_spec[1],3))
    Chz           = np.zeros((grid_spec[0],grid_spec[1],8))
    get_field_coefficients(Cex,Cey,Chz,grid_spec,params)
    
    Ex            = np.zeros((grid_spec[0],grid_spec[1]  ))
    Ey            = np.zeros((grid_spec[0],grid_spec[1]  ))
    Hz            = np.zeros((grid_spec[0],grid_spec[1],3))
    
    pulse         = get_pulse(grid_spec,params)
    
    print("Begin time updates ...")
    for n in range(params[0][1]):
        
        i,j       = pulse['loc']
        Hz[i,j,0] = Hz[i,j,0] + Chz[i,j,3]*pulse['val'][n]  # MAGNETIC SOURCE DISCONTINUITY
        
        for i,j in indeces['hz_inside' ]:                   # FOR PROPAGATION PROBLEMS, CORRESPONDS TO THE NON-PML REGION
            Hz[i,j,0] = Chz[i,j,0]*Hz[i,j,0] + Chz[i,j,1]*(Ex[i+1,j  ] - Ex[i,j]) + Chz[i,j,2]*(Ey[i,j+1] - Ey[i,j])
        for i,j in indeces['hz_outside']:                   # FOR PROPAGATION PROBLEMS, CORRESPONDS TO THE     PML REGION
            Hz[i,j,1] = Chz[i,j,4]*Hz[i,j,1] + Chz[i,j,5]*(Ey[i  ,j+1] - Ey[i,j]) 
            Hz[i,j,2] = Chz[i,j,6]*Hz[i,j,2] + Chz[i,j,7]*(Ex[i+1,j  ] - Ex[i,j])
            Hz[i,j,0] = Hz [i,j,1] + Hz[i,j,2]

        for i,j in indeces['ex']:
            Ex[i,j]   = Cex[i,j,0]*Ex[i,j] + Cex[i,j,1]*(Hz[i,j,0] - Hz[i-1,j  ,0])
 
        for i,j in indeces['ey']:
            Ey[i,j]   = Cey[i,j,0]*Ey[i,j] + Cey[i,j,1]*(Hz[i,j,0] - Hz[i  ,j-1,0])

        exportarray.append(extract_view_and_calculate_output(Ex,Ey,Hz[:,:,0],grid_spec,params))

    visualizer(exportarray,params,paths)
    print("Finished!")

if __name__ == "__main__":
    main()
