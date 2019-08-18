import h5py


def struct2dict(filename):
    with h5py.File(filename, 'r') as f:
        if len(f.keys()) > 2:
            print('Careful... you have more than one structure in this file')
            return
        struct_name = list(f.keys())[1]
        struct_dict = {}
        for key in f[struct_name].keys():
            struct_dict[key] = f[struct_name + '/' + key][:]
    return struct_dict
