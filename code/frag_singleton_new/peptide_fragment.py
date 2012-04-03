import sys, pybindgen
from pybindgen import *

mod = Module('PeptideFragmentSingleton')
mod.add_include('"peptide_fragment.hpp"')
mod.add_container('ListD', ReturnValue.new('double'), 'list')
c = mod.add_class('PeptideFragment')
c.add_constructor([
        param('int',    'massType', default_value="1"),
        ])
c.add_instance_attribute('a', 'ListD')
c.add_instance_attribute('b', 'ListD')
c.add_instance_attribute('c', 'ListD')
c.add_instance_attribute('x', 'ListD')
c.add_instance_attribute('y', 'ListD')
c.add_instance_attribute('z', 'ListD')
c.add_instance_attribute('zdot', 'ListD')
c.add_instance_attribute('mass', 'ListD')
#c.add_instance_attribute('pI', 'double')
#c.add_instance_attribute('composition', 'std::string')
c.add_method('composition',  retval('std::string'),  [])
c.add_method('sequence',     retval('std::string'),  [], is_const=True)
c.add_method('a_ions',       ReturnValue.new('ListD'),  [], is_virtual=True)
c.add_method('b_ions',       ReturnValue.new('ListD'),  [], is_virtual=True)
c.add_method('y_ions',       ReturnValue.new('ListD'),  [], is_virtual=True)
c.add_method('peptide_mass', ReturnValue.new('ListD'),  [], is_virtual=True)
c.add_method('pI',           ReturnValue.new('double'), [], is_virtual=True)
c.add_method('analyze', None, [
        param('char *', 'peptide'),
        param('int',    'charge',   default_value="1"),
        param('char *', 'modification',   default_value='(char *)""')
        ])
mod.generate(sys.stdout)
