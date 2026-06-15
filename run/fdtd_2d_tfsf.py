###################################################
### FINITE DIFFERENCE TIME DOMAIN (FDTD) 2D TEz ###
### SOLVER WITH:                                ###
###     TOTAL-FIELD/SCATTERED-FIELD EXCITATION  ###
###     SPLIT-PML ABC                           ###
### by: OAMEED NOAKOASTEEN                      ###
###################################################

import os
import numpy as np

from fdtd_1d import get_pulse              as get_pulse_aux
from fdtd_1d import get_field_coefficients as get_field_coefficients_aux

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
    c0                      = 1/np.sqrt(params[2][1]*params[2][0])
    lx                      = 1
    ly                      = lx
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

def get_tfsf_indeces(grid_spec):
    size_factor = 0.32
    thickness_x = int(size_factor*grid_spec[0])
    thickness_y = int(size_factor*grid_spec[1])
    origin      = (thickness_x + 1, thickness_y)
    indeces_O1  = []
    indeces_O2  = []
    indeces_O3  = []
    indeces_O4  = []
    indeces_I1  = []
    indeces_I4  = []
    for     j in range(thickness_x + 1, grid_spec[0] - thickness_x - 1):
        for i in [thickness_y - 1]:
            indeces_O1.append((i,j))
    for     j in [grid_spec[0] - thickness_x]:
        for i in range(thickness_y, grid_spec[1] - thickness_y):
            indeces_O2.append((i,j))
    for     j in range(thickness_x + 1, grid_spec[0] - thickness_x - 1):
        for i in [grid_spec[1] - thickness_y + 1]:
            indeces_O3.append((i,j))
    for     j in [thickness_x]:
        for i in range(thickness_y, grid_spec[1] - thickness_y):
            indeces_O4.append((i,j))
    for     j in range(thickness_x + 1, grid_spec[0] - thickness_x - 1):
        for i in [thickness_y]:
            indeces_I1.append((i,j))
    for     j in [thickness_x + 1]:
        for i in range(thickness_y, grid_spec[1] - thickness_y):
            indeces_I4.append((i,j))
    return {'O1':indeces_O1,'O2':indeces_O2,'O3':indeces_O3,'O4':indeces_O4,'I1':indeces_I1,'I4':indeces_I4, 'origin':origin}

def get_objects(MPex,MPey,MPhz,grid_spec,params):
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
    options     = {'noobjects'   : get_objects_noobjects   ,
                   'rectanglePEC': get_objects_rectanglePEC,
                   'rectangleDE' : get_objects_rectangleDE ,
                   'circlePEC'   : get_objects_circlePEC   ,
                   'circleDE'    : get_objects_circleDE     }    
    object_func = options[params[1]['name']]
    object_func(MPex,MPey,MPhz,grid_spec,params)

def get_pml_region(MPex,MPey,MPhz,grid_spec,params):
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

def get_grid_spec_aux(grid_spec):
    size_factor_pml       = 0.25
    size_factor_pulse_loc = 0.32
    length_factor         = 1 - size_factor_pml - size_factor_pulse_loc 
    thickness_x           = int(size_factor_pulse_loc*grid_spec[0])
    thickness_y           = int(size_factor_pulse_loc*grid_spec[1])
    start_x               = thickness_x + 1
    end_x                 = grid_spec[0] - start_x
    start_y               = thickness_y + 1
    end_y                 = grid_spec[0] - start_y
    length_x              = end_x - start_x
    length_y              = end_y - start_y
    length_d              = np.sqrt(np.power(length_x,2) + np.power(length_y,2))
    N                     = int(1.5*length_d/length_factor)
    return [N, grid_spec[2], grid_spec[3]]

def get_tfsf_connection_to_1d_aux(indeces_tfst,grid_spec,params):
    aux_size_factor      = 0.32
    aux_excitation_index = int(aux_size_factor*grid_spec[0])
    aux_slope            = params[0][3]*np.pi/180
    aux_unit_vector      = (1/np.sqrt(1 + np.power(aux_slope,2)), aux_slope/np.sqrt(1 + np.power(aux_slope,2)))
    bnd_origin           = indeces_tfst['origin']
    O1                   = []
    O2                   = []
    O3                   = []
    O4                   = []
    I1                   = []
    I4                   = []
    for i,j in indeces_tfst['O1']:
        vector           = (i - bnd_origin[0], (j + 1) - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        O1.append(aux_excitation_index + int(np.floor(projection)))
    O1 = [x for x in zip(indeces_tfst['O1'],O1)]
    for i,j in indeces_tfst['O2']:
        vector           = (i - bnd_origin[0],  j - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        O2.append(aux_excitation_index + int(np.floor(projection)))
    O2 = [x for x in zip(indeces_tfst['O2'],O2)]
    for i,j in indeces_tfst['O3']:
        vector           = (i - bnd_origin[0],  j - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        O3.append(aux_excitation_index + int(np.floor(projection)))
    O3 = [x for x in zip(indeces_tfst['O3'],O3)]
    for i,j in indeces_tfst['O4']:
        vector           = ((i + 1) - bnd_origin[0],  j - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        O4.append(aux_excitation_index + int(np.floor(projection)))
    O4 = [x for x in zip(indeces_tfst['O4'],O4)]
    for i,j in indeces_tfst['I1']:
        vector           = (i - bnd_origin[0],  (j - 1) - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        I1.append(aux_excitation_index + int(np.floor(projection)))
    I1 = [x for x in zip(indeces_tfst['I1'],I1)]
    for i,j in indeces_tfst['I4']:
        vector           = ((i - 1) - bnd_origin[0],  j - bnd_origin[1])
        projection       = np.dot(vector, aux_unit_vector)
        I4.append(aux_excitation_index + int(np.floor(projection)))
    I4 = [x for x in zip(indeces_tfst['I4'],I4)]
    return {'O1':O1,'O2':O2,'O3':O3,'O4':O4,'I1':I1,'I4':I4}

def get_field_discontinuity(ey,hz,index,region,params):
    def adjust_e_O1(ey,params):
        return ey*(-np.sin(params[0][3]*np.pi/180))
    def adjust_e_O2(ey,params):
        return ey*( np.cos(params[0][3]*np.pi/180))
    def adjust_e_O3(ey,params):
        return ey*(-np.sin(params[0][3]*np.pi/180))
    def adjust_e_O4(ey,params):
        return ey*( np.cos(params[0][3]*np.pi/180))
    def adjust_e_I1(ey,params):
        return ey
    def adjust_e_I4(ey,params):
        return ey
    options_adjust_e   = {'O1': adjust_e_O1,
                          'O2': adjust_e_O2,
                          'O3': adjust_e_O3,
                          'O4': adjust_e_O4,
                          'I1': adjust_e_I1,
                          'I4': adjust_e_I4 }
    adjust_e_func      = options_adjust_e[region]
    ey_value           = ey[index]
    hz_value           = hz[index]
    ey_value           = adjust_e_func(ey_value,params)
    return ey_value, hz_value

def extract_view_and_calculate_output(ex,ey,hz,grid_spec,params):
    def extract_view(ex,ey,hz,grid_spec):
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
    ex_copy, ey_copy, hz_copy = extract_view(ex,ey,hz,grid_spec)
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
    size_base  = 4
    cmap       = 'jet'
    fps        = 10
    dpi        = 90
    filename   = os.path.join(paths[0],'data' + '.gif')
    vmax       = get_color_maximum(exportarray,params) 
    fig, ax    = plt.subplots(figsize=(size_base,size_base))
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
    # params[0][3]: angle (deg) of wavefront
    # params[0][4]: problem type
    # params[1][*]: dictionay containing structure information
    # params[2][0]: FREE SPACE PERMITTIVITY
    # params[2][1]: FREE SPACE PERMEABILITY
    # params[2][2]: CFL FACTOR
    import argparse
    import shutil
    parser     = argparse.ArgumentParser()
    parser.add_argument('-p', type = str  , required = True)
    parser.add_argument('-v', type = str  , required = True)
    parser.add_argument('-n', type = int  , default  = 40  )
    parser.add_argument('-t', type = int  , default  = 500 )
    parser.add_argument('-l', type = int  , default  = 10  )
    parser.add_argument('-a', type = float, default  = 45  )
    args       = parser.parse_args()
    structure  = rJSON(args.p + '.json')
    params     = [[args.n,args.t,args.l,args.a,args.p],
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
    indeces_tfst  = get_tfsf_indeces  (grid_spec)
    
    Cex           = np.zeros((grid_spec[0],grid_spec[1],3))
    Cey           = np.zeros((grid_spec[0],grid_spec[1],3))
    Chz           = np.zeros((grid_spec[0],grid_spec[1],8))
    get_field_coefficients(Cex,Cey,Chz,grid_spec,params)
    
    Ex            = np.zeros((grid_spec[0],grid_spec[1]  ))
    Ey            = np.zeros((grid_spec[0],grid_spec[1]  ))
    Hz            = np.zeros((grid_spec[0],grid_spec[1],3))
    
    grid_spec_aux = get_grid_spec_aux(grid_spec)
    Cey_aux       = np.zeros((grid_spec_aux[0],3))
    Chz_aux       = np.zeros((grid_spec_aux[0],3))
    get_field_coefficients_aux(Cey_aux,Chz_aux,grid_spec_aux,params)
    
    Ey_aux        = np.zeros((grid_spec_aux[0],))
    Hz_aux        = np.zeros((grid_spec_aux[0],))
    
    indeces_tfst  = get_tfsf_connection_to_1d_aux(indeces_tfst,grid_spec_aux,params)
    
    pulse         = get_pulse_aux(grid_spec_aux,params)
    
    print("Begin time updates ...")
    for n in range(params[0][1]):
        
        # TFSF EXCITATION:
        # UPDATE OF THE 1D AUXILIARY GRID
        i         = pulse['loc']
        Hz_aux[i] = Hz_aux[i] + Chz_aux[i,2]*pulse['val'][n]
        
        for i in range(  grid_spec_aux[0]-1):
            Hz_aux[i] = Chz_aux[i,0]*Hz_aux[i] + Chz_aux[i,1]*(Ey_aux[i+1] - Ey_aux[i  ])

        for i in range(1,grid_spec_aux[0]-1):
            Ey_aux[i] = Cey_aux[i,0]*Ey_aux[i] + Cey_aux[i,1]*(Hz_aux[i  ] - Hz_aux[i-1])
        
        # TFSF EXCITATION:
        # ENFORCE CONSISTENCY UPDATES ALONG THE TFSF BOUNDARY
        for indeces_bnd, index_aux in indeces_tfst['O1']:
            e_disc, _      = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'O1',params)
            i,j            = indeces_bnd
            Hz[i,j,0]      = Hz[i,j,0] - Chz[i,j,1]*e_disc
        for indeces_bnd, index_aux in indeces_tfst['O2']:
            e_disc, h_disc = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'O2',params)
            i,j            = indeces_bnd
            Ey[i,j  ]      = Ey[i,j  ] + Cey[i,j,1]*h_disc
            Hz[i,j,0]      = Hz[i,j,0] + Chz[i,j,2]*e_disc
        for indeces_bnd, index_aux in indeces_tfst['O3']:
            e_disc, h_disc = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'O3',params)
            i,j            = indeces_bnd
            Ex[i,j  ]      = Ex[i,j  ] + Cex[i,j,1]*h_disc
            Hz[i,j,0]      = Hz[i,j,0] + Chz[i,j,1]*e_disc
        for indeces_bnd, index_aux in indeces_tfst['O4']:
            e_disc, _      = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'O4',params)
            i,j            = indeces_bnd
            Hz[i,j,0]      = Hz[i,j,0] - Chz[i,j,2]*e_disc
        for indeces_bnd, index_aux in indeces_tfst['I1']:
            _     , h_disc = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'I1',params)
            i,j            = indeces_bnd
            Ex[i,j  ]      = Ex[i,j  ] - Cex[i,j,1]*h_disc
        for indeces_bnd, index_aux in indeces_tfst['I4']:
            _     , h_disc = get_field_discontinuity(Ey_aux,Hz_aux,index_aux,'I4',params)
            i,j            = indeces_bnd
            Ey[i,j  ]      = Ey[i,j  ] - Cey[i,j,1]*h_disc
        
        # UPDATE OF THE 2D GRID
        for i,j in indeces['hz_inside' ]: # NON-PML REGION
            Hz[i,j,0] = Chz[i,j,0]*Hz[i,j,0] + Chz[i,j,1]*(Ex[i+1,j  ] - Ex[i,j]) + Chz[i,j,2]*(Ey[i,j+1] - Ey[i,j])
        for i,j in indeces['hz_outside']: #     PML REGION
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
