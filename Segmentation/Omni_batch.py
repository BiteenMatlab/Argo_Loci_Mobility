Phase_plot = False
Mask_plot = False

import numpy as np
from cellpose_omni import models, core

# for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 300
plt.style.use('dark_background')
# %matplotlib inline

# check if GPU is set up properly for segmentation 
use_GPU = core.use_gpu()
# print('>>> GPU activated? {}'.format(use_GPU))

# import Mask file list for segmentation 
file_list = open('mask_path_list.txt','r')
Mask_list = file_list.readlines()
Mask_list = [x.strip() for x in Mask_list]
# print(Mask_list)
# print(Mask_list[0])

from cellpose_omni import io, transforms
from omnipose.utils import normalize99
imgs = [io.imread(f) for f in Mask_list]
nimg = len(imgs)

if Phase_plot==True:
    # print some info about the images.
    for i in imgs:
        print('Original image shape:',i.shape)
        print('data type:',i.dtype)
        print('data range: min {}, max {}\n'.format(i.min(),i.max()))
    nimg = len(imgs)
    print('\nnumber of images:',nimg)

    fig = plt.figure(figsize=[40]*2,frameon=False) # initialize figure
    print('\n')
    for k in range(len(imgs)):
        img = transforms.move_min_dim(imgs[k]) # move the channel dimension last
        if len(img.shape)>2:
            # imgs[k] = img[:,:,1] # could pick out a specific channel
            imgs[k] = np.mean(img,axis=-1) # or just turn into grayscale 
            
        imgs[k] = normalize99(imgs[k])
        # imgs[k] = np.pad(imgs[k],10,'edge')
        print('new shape: ', imgs[k].shape)
        plt.subplot(1,len(Mask_list),k+1)
        plt.imshow(imgs[k],cmap='gray')
        plt.axis('off')


from cellpose_omni import models
from cellpose_omni.models import MODEL_NAMES
model_name = 'bact_phase_omni'
model = models.CellposeModel(gpu=use_GPU, model_type=model_name)

import time
chans = [0,0] #this means segment based on first channel, no second channel 

# n = [-1] # make a list of integers to select which images you want to segment
n = range(nimg) # or just segment them all 

# define parameters
params = {'channels':chans, # always define this with the model
          'rescale': None, # upscale or downscale your images, None = no rescaling 
          'mask_threshold': -1, # erode or dilate masks with higher or lower values 
          'flow_threshold': 0, # default is .4, but only needed if there are spurious masks to clean up; slows down output
          'transparency': True, # transparency in flow output
          'omni': True, # we can turn off Omnipose mask reconstruction, not advised 
          'cluster': True, # use DBSCAN clustering
          'resample': True, # whether or not to run dynamics on rescaled grid or original grid 
          # 'verbose': False, # turn on if you want to see more output 
          'tile': False, # average the outputs from flipped (augmented) images; slower, usually not needed 
          'niter': None, # None lets Omnipose calculate # of Euler iterations (usually <20) but you can tune it for over/under segmentation 
          'augment': False, # Can optionally rotate the image and average outputs, usually not needed 
          'affinity_seg': False, # new feature, stay tuned...
         }

tic = time.time() 
masks, flows, styles = model.eval([imgs[i] for i in n],**params)

net_time = time.time() - tic
print('total segmentation time: {}s'.format(net_time))

from cellpose_omni import plot
import omnipose

if Mask_plot==True:
    for idx,i in enumerate(n):

        maski = masks[idx] # get masks
        bdi = flows[idx][-1] # get boundaries
        flowi = flows[idx][0] # get RGB flows 

        # set up the output figure to better match the resolution of the images 
        f = 2
        szX = maski.shape[-1]/mpl.rcParams['figure.dpi']*f
        szY = maski.shape[-2]/mpl.rcParams['figure.dpi']*f
        fig = plt.figure(figsize=(szY,szX*4))
        fig.patch.set_facecolor([0]*4)
        
        plot.show_segmentation(fig, omnipose.utils.normalize99(imgs[i]), 
                            maski, flowi, bdi, channels=chans, omni=True)

        plt.tight_layout()
        plt.show()


io.save_masks(imgs, masks, flows, Mask_list, 
              png=True,
              tif=False, # whether to use PNG or TIF format
              suffix='', # suffix to add to files if needed 
              save_flows=False, # saves both RGB depiction as *_flows.png and the raw components as *_dP.tif
              save_outlines=False, # save outline images 
              dir_above=0, # save output in the image directory or in the directory above (at the level of the image directory)
              in_folders=True, # save output in folders (recommended)
              save_txt=False, # txt file for outlines in imageJ
              save_ncolor=False) # save ncolor version of masks for visualization and editing 