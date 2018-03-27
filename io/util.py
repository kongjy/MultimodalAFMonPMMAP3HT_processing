from igor.binarywave import load
from scipy.signal import detrend
import io
import numpy as np
import os

def load_ibw(path, flatten=True):
    """
    Given a path to an Igor Binary Wave, return the image file as a 3 dimensional 
    numpy array.
    
    Input:
        path: string file path to .ibw file
        flatten (optional): boolean input to flatten topography data.
    Output:
        data: 3 dimensional numpy array containing .ibw data.
    """
    data = load(path)['wave']['wData']

    # Flatten the topography data by extracting any linear response.
    if flatten == True:   
        flat_topo = data.copy()
        flat_topo[:, :, 0] = detrend(flat_topo[:, :, 0])
        data = flat_topo
        
    data = np.rot90(data)    
    return data

def hyperslice(hyper, start, stop, rows = None , cols = None):
    """
    Sums a range of wavenumbers within a hyperspectral image. 
    
    Paramters 
    ----------
        hyper: class from HyperImage
        start: int
            start wavenumber
        stop: int
            stop wavenumber
        rows: tuple, (start, stop)
    	    starting and ending indices for rows within image 
            be displayed. If all rows are desired, can 
            leave blank.  
        cols: tuple, (start, stop)
            same as above, but for columns. 
    Returns
    -------
        slc: ndarray
            sum of intensities between the specified start and
            stop wavenumbers 
        
    """
    #show entire hyper image if no tuples are passed into arguments
    if rows == None: 
        rows = (0,hyper.channel_data.shape[0])
    if cols == None:
        cols = (0, hyper.channel_data.shape[1])        

    wavenumber = hyper.wavelength_data.tolist()
    wavenumberlist = [int(x) for x in wavenumber]
    #flip start and stop indices because of the 
    #way the wavenumber data is stored. 
    start_index = wavenumberlist.index(stop)
    stop_index = wavenumberlist.index(start) 
    span = stop_index - start_index
    slc = hyper.hyper_image[rows[0]:rows[1],cols[0]:cols[1],start_index]
    for i in range(span):
        slc += hyper.hyper_image[rows[0]:rows[1],cols[0]:cols[1],start_index+i]
    
    return slc
    
    

class HyperImage():
    """
    A class representing a Hyper image. Give the path to the Hyper data,
    and receive a class that stores this information as a hyper image,
    and series of channel images.
    """
    def __init__(self, path):
        
        self.wavelength_data = None
        self.channel_names = []
        full_path = os.path.realpath(path)
        directory = os.path.dirname(full_path)
        
        # Get the scan parameters and channel details.
        self.parms, channels, e =  read_anfatec_params(full_path)
        
        x_pixel = int(self.parms['xPixel'])
        y_pixel = int(self.parms['yPixel'])
        
        self.wavelength_data = np.loadtxt(os.path.join(directory,str(channels[0]['FileNameWavelengths'])))
        wavenumber_length = self.wavelength_data.shape[0] 
        image_shape = (x_pixel,y_pixel,wavenumber_length)

        hyper_image = np.zeros(image_shape)
        
        # This scales the integer data into floats.
        pifm_scaling = float(channels[0]['Scale'])
        
        # Read the Raw Hyper data from the bitfile.
        data = np.fromfile(os.path.join(directory,channels[0]['FileName']),dtype='i4')
        for i,line in enumerate(np.split(data,y_pixel)):
            for j, pixel in enumerate(np.split(line,x_pixel)):
                    hyper_image[j,i,:] = pifm_scaling*pixel
                    
        # Put all the different channels into one big array.
        channel_data = np.zeros((x_pixel, y_pixel, len(channels[1:])))
        for ch, channel in enumerate(channels[1:]):
            self.channel_names.append(channel['Caption'])
            data = np.fromfile(os.path.join(directory,channel['FileName']),dtype='i4')
            scaling = float(channel['Scale'])

            for i,line in enumerate(np.split(data,y_pixel)):
                for j, pixel in enumerate(np.split(line,x_pixel)):
                        channel_data[j,i,ch] = (scaling*pixel)

        # Here's how we access the different hyper and channel data.
        self.hyper_image = np.rot90(hyper_image, k=-1)
        self.channel_data = np.rot90(channel_data, k=-1)
        
        self.hyper_image = self.hyper_image[:,::-1,:]
        self.channel_data = self.channel_data[:,::-1,:]
        

def read_anfatec_params(path):
    """
    Reads in an ANFATEC parameter file. This file is produced by the Molecular
    Vista PiFM system and describes all parameters need to interpret the data 
    files produced when the data is saved.
    
    Input:
        path: a path to the ANFATEC parameter file.
        
    Output:
        file_descriptions: A list of dictionaries, with each item in the list 
            corresponding to a channel that was recorded by the PiFM.
        scan_params: A dictionary of non-channel specific scan parameters.
        
    """
    file_descriptions = []
    spectra_descriptions = []
    scan_params = {}
    parameters = {}
    inside_description = False


    with io.open(path,  'r', encoding = "ISO-8859-1") as f:
        
        for i,row in enumerate(f): 
            
            # Get rid of newline characters at the end of the line.
            row = row.strip()
            #check to make sure its  not empty 
            if row:
                # First line of the file is useless. We tell the reader to stop at ';'
                if row[0] == unicode(';'):
                    continue
                
                # This string indicates that we have reached a channel description.
                if row.endswith('Begin') & row.startswith('File'):
                    inside_description = True
                    continue
                if row.endswith('SpectrumDescBegin'):
                    inside_description = True
                    continue
                if row.endswith('End') & row.startswith('File'):
                    file_descriptions.append(parameters)
                    parameters = {}
                    inside_description = False
                if row.endswith('SpectrumDescEnd'):
                    spectra_descriptions.append(parameters)
                    parameters = {}
 
                #split between :; creates list of two elements 
                split_row = row.split(':')
                
                for i, el in enumerate(split_row):
                    split_row[i] = el.strip()
                
                # We want to save the channel parameters to a separate structure.
                if inside_description:
                    parameters[split_row[0]] = split_row[-1]
                else:
                    scan_params[split_row[0]] = split_row[-1]
            
                    
    return scan_params, file_descriptions, spectra_descriptions

def to_2d(spectralmatrix, laserpower):
    """
    Takes a raw 3-D N x M x S matrix and restructures it into a laser power corrected,
    2D (NxM) x S matrix. 
    """
    xaxis=spectralmatrix.shape[0]
    yaxis=spectralmatrix.shape[1]
    zaxis=spectralmatrix[0,0,:].shape[0]
    featurematrix=np.zeros((xaxis*yaxis, zaxis))
    counter=-1
    #loop over elements in spectralmatrix 
    for x in range(xaxis):
        for y in range(yaxis):
            counter+=1 
            for z in range(zaxis):
                featurematrix[counter]=(spectralmatrix[x,y,:])/laserpower
    return featurematrix