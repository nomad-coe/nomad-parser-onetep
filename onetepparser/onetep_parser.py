#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from nomad.units import ureg
from nomad.parsing.file_parser import BasicParser


class OnetepParser(BasicParser):
    def __init__(self):
        re_f = r'\-*\d+\.\d+E*e*\-*\+*\d*'

        super().__init__(
            specifications=dict(
                name='parsers/onetep', code_name='ONETEP', code_homepage='https://www.onetep.org/',
                domain='dft', mainfile_contents_re=r'####### #     # ####### ####### ####### ######'),
            units_mapping=dict(energy=ureg.hartree, length=ureg.bohr),
            auxilliary_files=r'([\w\-]+\.dat)',
            program_version=r'Version\s*([\d\.]+)',
            lattice_vectors=r'\%block lattice_cart\s*([\s\S]+?)\%endblock lattice_cart',
            atom_labels_atom_positions=rf'\%block positions\_abs\s*(\w+\s+{re_f}\s+{re_f}\s+{re_f}[\s\S]+?)\%endblock positions\_abs',
            XC_functional=r'xc\_functional\s*\:\s*(\w+)',
            energy_total=rf'Total energy\s*=\s*({re_f})\s*Eh',
            atom_forces=r'Forces \*+([\s\S]+?)\*{50}')
