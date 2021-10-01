# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from pwem.convert import AtomicStructHandler

def generate_gpf(protFile, spacing, xc, yc, zc, npts, outDir, ligandFns=None):
  """
  Build the GPF file that is needed for AUTOGRID to generate the electrostatic grid
  """

  protName, protExt = os.path.splitext(os.path.basename(protFile))
  gpf_file = os.path.join(outDir, protName + '.gpf')
  npts = round(npts, None)

  protAtomTypes = parseAtomTypes(protFile)
  protAtomTypes = ' '.join(sortSet(protAtomTypes))

  if ligandFns == None:
      ligAtomTypes = 'A C HD N NA OA SA'
  else:
      ligAtomTypes = set([])
      for ligFn in ligandFns:
          ligAtomTypes = ligAtomTypes | parseAtomTypes(ligFn)

      ligAtomTypes = ' '.join(sortSet(ligAtomTypes))



  with open(os.path.abspath(gpf_file), "w") as file:
    file.write("npts %s %s %s                        # num.grid points in xyz\n" % (npts, npts, npts))
    file.write("gridfld %s.maps.fld                # grid_data_file\n" % (protName))
    file.write("spacing %s                          # spacing(A)\n" % (spacing))
    file.write("receptor_types %s     # receptor atom types\n" % (protAtomTypes))
    file.write("ligand_types %s       # ligand atom types\n" % (ligAtomTypes))
    file.write("receptor %s                  # macromolecule\n" % (os.path.abspath(protFile)))
    file.write("gridcenter %s %s %s           # xyz-coordinates or auto\n" % (xc, yc, zc))
    file.write("smooth 0.5                           # store minimum energy w/in rad(A)\n")
    for ligType in ligAtomTypes.split():
        file.write("map %s.%s.map                       # atom-specific affinity map\n" % (protName, ligType))
    file.write("elecmap %s.e.map                   # electrostatic potential map\n" % (protName))
    file.write("dsolvmap %s.d.map                  # desolvation potential map\n" % (protName))
    file.write("dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant")

  return os.path.abspath(gpf_file)

def parseAtomTypes(pdbFile):
    atomTypes = set([])
    with open(pdbFile) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
              pLine = line.split()
              at = pLine[12]
              atomTypes.add(at)
    return atomTypes

def calculate_centerMass(atomStructFile):
      """
      Returns the geometric center of mass of an Entity (anything with a get_atoms function in biopython).
      Geometric assumes all masses are equal
      """

      try:
          structureHandler = AtomicStructHandler()
          structureHandler.read(atomStructFile)
          center_coord = structureHandler.centerOfMass(geometric=True)
          structure = structureHandler.getStructure()

          return structure, center_coord[0], center_coord[1], center_coord[2]  # structure, x,y,z

      except Exception as e:
          print("ERROR: ", "A pdb file was not entered in the Atomic structure field. Please enter it.", e)
          return
      
def sortSet(seti):
    seti = list(seti)
    seti.sort()
    return seti