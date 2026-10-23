import numpy as np
import ctypes as ct
import tables as tb
import os

# We have to do a dirty trick to get the correct path to the native upside library.
# It is likely possible to modify the build system to drop a config file 
# that contains the necessary paths.  Since I do a lot of work concurrently
# in the source tree, I don't want to depend on the build system to copy the
# python executables into obj/.  This would make it easy not to commit changes
# back when live-editing.

# FIXME This could fail if the user links directly into the py directory
# instead of going through upside/py.  
py_source_dir = os.path.dirname(__file__)
obj_dir = os.path.join(py_source_dir, '..', 'obj')


def load_upside_library():
    candidate_names = [
        'libupside.so',
        'libupside.dylib',
        'libupside.dll',
    ]
    attempted_paths = []

    for library_name in candidate_names:
        library_path = os.path.join(obj_dir, library_name)
        attempted_paths.append(library_path)
        if os.path.exists(library_path):
            return ct.cdll.LoadLibrary(library_path)

    raise OSError(
        'Could not find the Upside native library. Tried: {}'.format(
            ', '.join(attempted_paths)
        )
    )


calc = load_upside_library()

calc.construct_deriv_engine.restype  = ct.c_void_p
calc.construct_deriv_engine.argtypes = [ct.c_int, ct.c_char_p]

calc.free_deriv_engine.restype  = None
calc.free_deriv_engine.argtypes = [ct.c_void_p]

calc.evaluate_energy.restype  = ct.c_int
calc.evaluate_energy.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p]

calc.evaluate_deriv.restype  = ct.c_int
calc.evaluate_deriv.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p]

calc.set_param.restype  = ct.c_int
calc.set_param.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_param_deriv.restype  = ct.c_int
calc.get_param_deriv.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_param.restype  = ct.c_int
calc.get_param.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_param_size.restype  = ct.c_int
calc.get_param_size.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_output_dims.restype  = ct.c_int
calc.get_output_dims.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_sens.restype  = ct.c_int
calc.get_sens.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_output.restype  = ct.c_int
calc.get_output.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p]

calc.get_value_by_name.restype  = ct.c_int
calc.get_value_by_name.argtypes = [ct.c_int, ct.c_void_p, ct.c_void_p, ct.c_char_p, ct.c_char_p]


def in_process_upside(args, verbose=True):
    # A bit of paranoia to make sure this doesn't interfere with python handlers
    # I more sophisticated strategy would modify the signal handlers to interoperate with
    # the Python handlers
    # Additionally, add a first argument because the program name is expected in argv[0]
    exec_args = ['python_library', '--re-raise-signal'] + args
    exec_c_strs = [ct.c_char_p(s.encode('utf-8')) for s in exec_args]

    # c_str_length = 1+max(len(a) for a in args)
    # print 'arg length %i' % c_str_length

    # build the strings into a char** for C
    c_arr_t = (ct.c_char_p * len(exec_args))
    c_arr = c_arr_t()
    for i,a in enumerate(exec_c_strs):
        c_arr[i] = a

    calc.upside_main.restype  = ct.c_int
    # calc.upside_main.argtypes = [ct.c_int, ct.POINTER(ct.c_char_p)]
    calc.upside_main.argtypes = [ct.c_int, c_arr_t, ct.c_int]

    retcode = calc.upside_main(len(exec_args), c_arr, int(verbose))
    if retcode:
        raise RuntimeError('In process Upside returned %i'%retcode)
    return


class Upside(object):
    def __init__(self, config_file_path):
        self.config_file_path = str(config_file_path)
        with tb.open_file(self.config_file_path) as t:
            self.initial_pos = t.root.input.pos[:,:,0]
            self.n_atom = self.initial_pos.shape[0]
        self.engine = calc.construct_deriv_engine(self.n_atom, bytes(self.config_file_path, encoding="ascii"))
        if self.engine is None: raise RuntimeError('Unable to initialize upside engine for %s'%(config_file_path,))

    def __repr__(self):
        return 'Upside(%r, %r)'%(self.n_atom, self.config_file_path)

    def energy(self, pos):
        pos = np.require(pos, dtype='f4', requirements='C')
        assert pos.shape == (self.n_atom,3)
        energy = np.zeros(1, dtype='f4')
        retcode = calc.evaluate_energy(energy.ctypes.data, self.engine, pos.ctypes.data)
        if retcode: raise RuntimeError('Unable to evaluate energy')
        return energy[0]

    def deriv(self, pos):
        pos = np.require(pos, dtype='f4', requirements='C')
        assert pos.shape == (self.n_atom,3)
        deriv = np.zeros_like(pos)
        retcode = calc.evaluate_deriv(deriv.ctypes.data, self.engine, pos.ctypes.data)
        if retcode: raise RuntimeError('Unable to evaluate derivative')
        return deriv

    def set_param(self, param, node_name):
        node_name = bytes(node_name, encoding="ascii")
        param_size = param.shape
        param = np.require(param.ravel(), dtype='f4', requirements='C')  # flatten and make contiguous
        retcode = calc.set_param(int(param.shape[0]), param.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to set param with size %s for node %s'%(param_size,node_name))

    def get_param_deriv(self, param_shape, node_name):
        node_name = bytes(node_name, encoding="ascii")
        n_param = int(np.prod(param_shape))
        deriv = np.zeros(param_shape, dtype='f4')
        retcode = calc.get_param_deriv(n_param, deriv.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get param deriv')
        return deriv

    def get_param(self, param_shape, node_name):
        node_name = bytes(node_name, encoding="ascii")
        n_param = int(np.prod(param_shape))
        param = np.zeros(param_shape, dtype='f4')
        retcode = calc.get_param(n_param, param.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get param')
        return param

    def get_param_size(self, node_name):
        node_name = bytes(node_name, encoding="ascii")
        n_param = np.zeros(1, dtype=np.intc)
        retcode = calc.get_param_size(n_param.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get param size')
        return int(n_param[0])

    def get_sens(self, node_name):
        node_name = bytes(node_name, encoding="ascii")
        n_elem = np.zeros(1,dtype=np.intc)
        elem_width = np.zeros(1,dtype=np.intc)
        retcode = calc.get_output_dims(n_elem.ctypes.data, elem_width.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get output dims')
        output_shape = (int(n_elem[0]), int(elem_width[0]))

        n_output = int(np.prod(output_shape))
        output = np.zeros(output_shape, dtype='f4')
        retcode = calc.get_sens(n_output, output.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get output')
        return output

    def get_output(self, node_name):
        node_name = bytes(node_name, encoding="ascii")
        n_elem = np.zeros(1,dtype=np.intc)
        elem_width = np.zeros(1,dtype=np.intc)
        retcode = calc.get_output_dims(n_elem.ctypes.data, elem_width.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get output dims')
        output_shape = (int(n_elem[0]), int(elem_width[0]))

        n_output = int(np.prod(output_shape))
        output = np.zeros(output_shape, dtype='f4')
        retcode = calc.get_output(n_output, output.ctypes.data, self.engine, node_name)
        if retcode: raise RuntimeError('Unable to get output')
        return output

    def get_value_by_name(self, value_shape, node_name, log_name):
        node_name = bytes(node_name, encoding="ascii")
        log_name = bytes(log_name, encoding="ascii")
        n_param = int(np.prod(value_shape))
        value = np.zeros(value_shape, dtype='f4')
        retcode = calc.get_value_by_name(n_param, value.ctypes.data, self.engine, node_name, log_name)
        if retcode: raise RuntimeError('Unable to get value by name')
        return value

    def __del__(self):
        calc.free_deriv_engine(self.engine)
