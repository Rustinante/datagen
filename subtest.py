import h5py
import multiprocessing
from multiprocessing import Queue
import numpy as np


class HDF5FileLoaderSubprocess(multiprocessing.Process):
    def __init__(self, filename, feature_dict_key, label_dict_key, output_queue):
        """
        :type output_queue: Queue
        :param filename:
        :param feature_dict_key:
        :param label_dict_key:
        :param output_queue:
        """
        super(HDF5FileLoaderSubprocess, self).__init__()
        
        self.filename = filename
        self.feature_dict_key = feature_dict_key
        self.label_dict_key = label_dict_key
        self.output_queue = output_queue
        
        self.hdf5_file = h5py.File(filename, 'r')
        self.features = self.hdf5_file[feature_dict_key]
        self.labels = self.hdf5_file[label_dict_key]
        self.length = self.labels.shape[0]
        self.permutation_index = np.arange(self.length)
        
    def run(self):
        while True:
            np.random.shuffle(self.permutation_index)
            for index in self.permutation_index:
                data_point = (self.features[index], self.labels[index])
                self.output_queue.put(data_point)
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, item):
        return self.features[item], self.labels[item]
    
    def __iter__(self):
        for index in range(self.length):
            yield self.features[index], self.labels[index]
        
    def __del__(self):
        self.hdf5_file.close()


q=Queue(4)
s=HDF5FileLoaderSubprocess('chr20_train.align.hdf5', 'feature/data', 'label/data', q)