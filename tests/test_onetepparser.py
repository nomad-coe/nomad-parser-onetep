#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest

from nomad.datamodel import EntryArchive
from onetepparser import OnetepParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return OnetepParser()


def test_basic(parser):
    archive = EntryArchive()

    parser.parse('tests/data/test05/test05.out', archive, None)

    sec_run = archive.run[0]
    assert sec_run.program.version == '4.4.1'

    sec_system = archive.run[0].system[0]
    assert sec_system.atoms.lattice_vectors[2][2].magnitude == approx(2.0108734e-09)
    assert sec_system.atoms.labels[3] == 'Cl'
    assert sec_system.atoms.positions[4][2].magnitude == approx(1.01078267e-09)

    sec_scc = sec_run.calculation[0]
    assert sec_scc.energy.total.value.magnitude == approx(-2.24380887e-16)


def test_1(parser):
    archive = EntryArchive()

    parser.parse('tests/data/test08/ethene_relax.out', archive, None)

    sec_sccs = archive.run[0].calculation
    assert len(sec_sccs) == 33
    assert sec_sccs[8].energy.total.value.magnitude == approx(-5.95750082e-17)
    assert sec_sccs[2].forces.total.value[5][1].magnitude == approx(-6.68217323e-10)
