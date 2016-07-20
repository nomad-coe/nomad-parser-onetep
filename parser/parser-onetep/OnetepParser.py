from __future__ import division
from builtins import str
from builtins import range
from builtins import object
import setup_paths
import numpy as np
import math
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from OnetepCommon import get_metaInfo
import OnetepCellParser
import OnetepBandParser
import OnetepMDParser
import OnetepTSParser
import logging, os, re, sys


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  PARSER CONTEXT CLASS  ##################################################################################################################
################################################################################################################################################################
############################################################################ onetep.Parser Version 1.0 #########################################################
################################################################################################################################################################

logger = logging.getLogger("nomad.OnetepParser")

class OnetepParserContext(object):

    def __init__(self):
        """ Initialise variables used within the current superContext """
        self.functionals                       = []
        self.func_total                        = []
        self.relativistic                      = []
        self.dispersion                        = None
        self.cell                              = []
        self.at_nr                             = 0
        self.atom_type_mass                    = []
        self.atom_labels                       = []
        self.atom_optim_position = []
        self.onetep_optimised_atom_positions   = [] 
        self.onetep_atom_optim_label           = []
        self.atom_forces                       = []
        self.atom_forces_band                       = []
        self.stress_tensor_value               = []
        self.raman_tensor_value               = []
        self.onetep_atom_positions              = []
        self.atom_positions                     = []
        self.a                                 = []
        self.b                                 = []
        self.c                                 = []
        self.alpha                             = []
        self.beta                              = []
        self.gamma                             = []
        self.volume                            = 0
        self.time_calc                         = 0
        self.time_calculation                  = []
        self.energy_total_scf_iteration_list   = []
        self.scfIterNr                         = []
        self.functional_type                   = []
        self.functional_weight                 = None
        self.k_nr                              = 0
        self.e_nr                              = 0
        self.k_count_1                         = 0
        self.k_nr_1                            = 0
        self.e_nr_1                            = 0
        self.onetep_band_kpoints               = []
        self.onetep_band_energies              = []
        self.onetep_band_kpoints_1             = []
        self.onetep_band_energies_1            = []
        self.k_path_nr                         = 0
        self.band_en                           = []
        self.van_der_waals_name                = None
        self.e_spin_1                          = []
        self.e_spin_2                          = []
        self.optim_method                      = []
        self.disp_energy                       = [] 
        self.geoConvergence = None
        self.total_energy_corrected_for_finite_basis = []
        self.contr_s = []
        self.contr_p = []
        self.contr_d = []
        self.contr_f = []
        self.total_contribution = []
        self.total_charge = []
        self.frequencies = []
        self.ir_intens = []
        self.raman_act = []
        self.irr_repres = []
        self.nr_iter =[]
        self.onetep_atom_velocities = []
        self.time_0 = []
        self.frame_temp = []
        self.frame_press =[]
        self.frame_potential =[]
        self.frame_total_energy = []
        self.frame_hamiltonian = []
        self.frame_kinetic = []
        self.frame_atom_forces=[]
        self.frame_atom_lables =[]
        self.frame_stress_tensor =[]
        self.frame_position =[]
        self.frame_cell=[]
        self.energy_frame_gain = []
        self.frame_time=[]
        self.energy_frame =[]
        self.total_energy_frame = []
        self.frame_energies =[]
        self.scfgIndex=[]
        self.n_spin_channels = []
        self.n_spin_channels_bands = []
        self.energy_frame_free = []
        self.energy_frame_T0 = []
        self.scf_conv_thresh = []
        self.n_iteration = []
        self.wall_time_end =[]
        self.initial_scf_iter_time = []
        self.ts_total_energy = []
        self.ts_cell_vector = []
        self.ts_forces = []
        self.ts_positions = []
        self.ts_path =[]
        self.ts_total_energy_f = []
        self.ts_cell_vector_f = []
        self.ts_forces_f = []
        self.ts_positions_f = []
        self.ts_path_f =[]
        self.ts_total_energy_p = []
        self.ts_cell_vector_p = []
        self.ts_forces_p = []
        self.ts_positions_p = []
        self.ts_path_p =[]
        self.total_spin_mulliken = []
        self.kinetic = []
        self.local_pseudo = [] 
        self.nonlocal_pseudo = []
        
        self.Hartree = []
        self.Exchangecorrelation = [] 
        self.Ewald = []
        self.DispersionCorrection =[] 
        self.Integrateddensity = []
    def initialize_values(self):
        """ Initializes the values of variables in superContexts that are used to parse different files """
        self.pippo = None
        self.at_nr_opt = None
    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()
    
    def onClose_section_topology(self, backend, gIndex, section):
    # Processing the atom masses and type
        #get cached values of onetep_store_atom_mass and type
        mass = section['x_onetep_store_atom_mass']
        name = section['x_onetep_store_atom_name']
        
        n_ngwf = section['x_onetep_n_ngwf_atom_store']
        radius = section['x_onetep_ngwf_radius_store']
        if mass:
            for i in range(len(mass)):
                #converting mass in Kg from AMU
                # mass[i] = float(mass[i]) * 1.66053904e-27

                backend.openSection('section_atom_type')
                backend.addValue('atom_type_mass', float(mass[i]))
                backend.addValue('atom_type_name', name[i])
                backend.closeSection('section_atom_type', gIndex+i)
                backend.addValue('x_onetep_n_ngwf_atom', n_ngwf[i])
                backend.addValue('x_onetep_ngwf_radius', radius[i])
# Translating the XC functional name to the NOMAD standard
    def onClose_x_onetep_section_functionals(self, backend, gIndex, section):
        """When all the functional definitions have been gathered, matches them
        with the nomad correspondents and combines into one single string which
        is put into the backend.
        """
        # Get the list of functional and relativistic names
        functional_names = section["x_onetep_functional_name"]
        
        #relativistic_names = section["onetep_relativity_treatment_scf"]

        
        functional_map = {
            "PBE": "GGA_C_PBE_GGA_X_PBE",
            "LDA": "LDA_C_PZ_LDA_X_PZ",
            "PW91": "GGA_C_PW91_GGA_X_PW91",
            "PW92": "GGA_C_PW92_GGA_X_PW92",
            "GGA": "GGA_X_PBE",
            "CAPZ": "LDA_C_PZ",
            "VWN": "LDA_C_VWN",
            "PBESOL": "GGA_X_PBE_SOL",          
            "RPBE": "GGA_X_RPBE",
            "REVPBE":"GGA_X_PBE_R",
            "BLYP":"GGA_X_B88_GGA_C_LYP",
            "XLYP":"GGA_XC_XLYP",
            "OPTB88":"GGA_X_OPTB88_VDW_LDA_C_PZ",
            "OPTPBE":"GGA_X_OPTPBE_VDW_LDA_C_PZ",
            "VDWDF":"GGA_X_PBE_R_VDW_LDA_C_PZ",
            "VDWDF2":"GGA_X_RPW86_LDA_C_PZ",
            "VV10":"HYB_GGA_XC_LC_VV10",
            "AVV10S":"GGA_X_AM05_GGA_C_AM05_rVV10-sol",
            
        }
       
       
        # Match each onetep functional name and sort the matches into a list
        self.functionals = []

        for name in functional_names:
            match = functional_map.get(name)
            if match:
                self.functionals.append(match)
            
        # self.functionals = "_".join(sorted(self.functionals))
        

        
    def onClose_x_onetep_section_relativity_treatment(self, backend, gIndex, section):
        relativistic_names = section["x_onetep_relativity_treatment_scf"]
        
        # Define a mapping for the relativistic treatments
        relativistic_map = {
            " Koelling-Harmon": "scalar_relativistic"
        }
        self.relativistic = []

        for name in relativistic_names:
            match = relativistic_map.get(name)
            if match:
                self.relativistic.append(match)
        self.relativistic = "_".join(sorted(self.relativistic))
# Here we add info about the XC functional and relativistic treatment
    
    def onClose_x_onetep_section_geom_optimisation_method(self, backend, gIndex, section):
        self.optim_method = section["x_onetep_geometry_optim_method"]

    def onClose_section_method(self, backend, gIndex, section):
        
        self.van_der_waals_name = section["van_der_Waals_method"]
        
        if self.van_der_waals_name is not None:
            dispersion_map = {
            " OBS correction scheme": "OBS",
            " G06 correction scheme": "G06",
            " TS correction scheme" : "TS",
            " JCHS correction scheme" : "JCHS",
            }
            self.dispersion = []
            for name in self.van_der_waals_name:
                match = dispersion_map.get(name)
                if match:
                    self.dispersion.append(match)    
            self.dispersion = "_".join(sorted(self.dispersion))
        else:
            pass
        #backend.addValue('van_der_Waals_method',self.dispersion)
        
        
        if self.functional_weight is not None:
            self.func_and_weight = []
            for i in range(len(self.functional_types)):
                self.func_total.append(self.functionals[i]+'_'+self.functional_weight[i]) 
                backend.openSection('section_XC_functionals')            
                backend.addValue('XC_functional_name', self.functionals[i]) 
                backend.addValue('XC_functional_weight', self.functional_weight[i])
                backend.closeSection('section_XC_functionals',gIndex+i)        
        # Push the functional string into the backend
        # Push the relativistic treatment string into the backend
            backend.addValue('XC_functional', "_".join(sorted(self.functionals)))
            # backend.addValue('relativity_method', self.relativistic)
        #if self.functional_weight = 0
            if self.dispersion is not None:
                backend.addValue('XC_method_current', ("_".join(sorted(self.func_total)))+'_'+self.dispersion+'_'+self.relativistic)
            else:
                backend.addValue('XC_method_current', ("_".join(sorted(self.func_total)))+'_'+self.relativistic)
        else:
            for i in range(len(self.functionals)):
          #      self.func_total.append(self.functionals[i]+'_'+self.functional_weight[i])
                backend.openSection('section_XC_functionals')            
                backend.addValue('XC_functional_name', self.functionals[i]) 
         #       backend.addValue('XC_functional_weight', self.functional_weight[i])
                backend.closeSection('section_XC_functionals',gIndex+i)
            backend.addValue('XC_functional', "_".join(sorted(self.functionals)))
            # backend.addValue('relativity_method', self.relativistic)
            if self.dispersion is not None:
                backend.addValue('XC_method_current', ("_".join(sorted(self.functionals)))+'_'+self.dispersion+'_'+self.relativistic)
            else:
                backend.addValue('XC_method_current', ("_".join(sorted(self.functionals))))
        
        if self.n_spin_channels:
            backend.addValue('number_of_spin_channels', self.n_spin_channels[0])
         

    
    
# Here we add basis set name and kind for the plane wave code
    def onClose_section_basis_set_cell_dependent(self, backend, gIndex, section):
        ecut_str = section['x_onetep_basis_set_planewave_cutoff']
        self.ecut = float(ecut_str[0])
        eVtoRy = 0.073498618
        ecut_str_name = int(round(eVtoRy*self.ecut))

        basis_set_kind = 'psinc_functions'
        basis_set_name = 'psinc_functions'
        backend.addValue('basis_set_planewave_cutoff', self.ecut)
        backend.addValue('basis_set_cell_dependent_kind', basis_set_kind)
        backend.addValue('basis_set_cell_dependent_name', basis_set_name)

    def onClose_x_onetep_section_cell_optim(self, backend, gIndex, section):
        """trigger called when _onetep_section_cell is closed"""
        # get cached values for onetep_cell_vector
        vet = section['x_onetep_cell_vector_optim']
     
        vet[0] = vet[0].split()
        vet[0] = [float(i) for i in vet[0]]

        vet[1] = vet[1].split()
        vet[1] = [float(i) for i in vet[1]]

        vet[2] = vet[2].split()
        vet[2] = [float(i) for i in vet[2]]

        self.cell.append(vet[0])
        self.cell.append(vet[1])
        self.cell.append(vet[2]) # Reconstructing the unit cell vector by vector    

     
# # Here we recover the unit cell dimensions (both magnitudes and angles) (useful to convert fractional coordinates to cartesian)
#     def onClose_x_onetep_section_atom_positions(self, backend, gIndex, section):
#         """trigger called when _onetep_section_atom_positions is closed"""
#         # get cached values for cell magnitudes and angles
#         self.a = section['x_onetep_cell_length_a']
#         self.b = section['x_onetep_cell_length_b']
#         self.c = section['x_onetep_cell_length_c']
#         self.alpha = section['x_onetep_cell_angle_alpha']
#         self.beta  = section['x_onetep_cell_angle_beta']
#         self.gamma = section['x_onetep_cell_angle_gamma']
#         self.volume = np.sqrt( 1 - math.cos(np.deg2rad(self.alpha[0]))**2
#                                  - math.cos(np.deg2rad(self.beta[0]))**2
#                                  - math.cos(np.deg2rad(self.gamma[0]))**2
#                                  + 2 * math.cos(np.deg2rad(self.alpha[0]))
#                                      * math.cos(np.deg2rad(self.beta[0]))
#                                      * math.cos(np.deg2rad(self.gamma[0])) ) * self.a[0]*self.b[0]*self.c[0]

    def onClose_x_onetep_section_atom_positions_optim(self, backend, gIndex, section):
        """trigger called when _onetep_section_atom_positions is closed"""
        # get cached values for cell magnitudes and angles
    
        self.a = section['x_onetep_cell_length_a_optim']
        self.b = section['x_onetep_cell_length_b_optim']
        self.c = section['x_onetep_cell_length_c_optim']
        self.alpha = section['x_onetep_cell_angle_alpha_optim']
        self.beta  = section['x_onetep_cell_angle_beta_optim']
        self.gamma = section['x_onetep_cell_angle_gamma_optim']
        self.volume = np.sqrt( 1 - math.cos(np.deg2rad(self.alpha[0]))**2
                                 - math.cos(np.deg2rad(self.beta[0]))**2
                                 - math.cos(np.deg2rad(self.gamma[0]))**2
                                 + 2 * math.cos(np.deg2rad(self.alpha[0]))
                                     * math.cos(np.deg2rad(self.beta[0]))
                                     * math.cos(np.deg2rad(self.gamma[0])) ) * self.a[0]*self.b[0]*self.c[0]    
       
# Storing the total energy of each SCF iteration in an array
    

# Processing forces acting on atoms (final converged forces)
    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        self.time_0 = section['x_onetep_frame_time_0'] 
        
        self.initial_scf_iter_time = section['x_onetep_initial_scf_iteration_wall_time']
        
        f_st = section['x_onetep_store_atom_forces']
        
        if f_st is not None:# and self.at_nr_opt is not None:
            for i in range(0, len(f_st)):
                f_st[i] = f_st[i].split()
                f_st[i] = [float(j) for j in f_st[i]]
                f_st_int = f_st[i]
                     
                self.atom_forces.append(f_st_int)
                    
                self.atom_forces = self.atom_forces[-len(f_st):] 
                    
            backend.addArrayValues('atom_forces', np.asarray(self.atom_forces))
        else:   
        #     for i in range(0, self.at_nr_opt):
        #         f_st[i] = f_st[i].split()
        #         f_st[i] = [float(j) for j in f_st[i]]
        #         f_st_int = f_st[i]
                     
        #         self.atom_forces.append(f_st_int)
                    
        #         self.atom_forces = self.atom_forces[-self.at_nr_opt:] 
                    
        #     backend.addArrayValues('atom_forces', np.asarray(self.atom_forces))
        # elif f_st is not None:
            
            pass        
# Add SCF k points and eigenvalue from *.band file to the backend (ONLY FOR SINGLE POINT CALCULATIONS AT THIS STAGE)
        for i in range(len(self.n_spin_channels)):
            if self.n_spin_channels[i]==1:
  
                backend.openSection('section_eigenvalues') # opening first section_eigenvalues
                backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
                backend.addArrayValues('eigenvalues_values', np.asarray(self.e_spin_1))
                backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
                backend.addValue('number_of_eigenvalues', self.e_nr_scf)
                backend.closeSection('section_eigenvalues', gIndex)
       
            if self.n_spin_channels[i]==2:         
                # print self.n_spin_channels,'ciao'
                # self.e_spin_2.(self.e_spin_1)
                both_spin = [[self.e_spin_2],[self.e_spin_1]]
                # 
                backend.openSection('section_eigenvalues') # opening the second section_eigenvalues (only for spin polarised calculations)
                backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
                backend.addArrayValues('eigenvalues_values', np.asarray(both_spin))
                backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
                backend.addValue('number_of_eigenvalues', self.e_nr_scf)

                backend.closeSection('section_eigenvalues', gIndex)    
        
        
        #backend.openSection('section_stress_tensor')
        # if self.stress_tensor_value:
            # backend.addArrayValues('stress_tensor',np.asarray(self.stress_tensor_value))
        #backend.closeSection('section_stress_tensor', gIndex)
        #backend.addValue('time_calculation', self.time_calc)

        #######Total Dispersion energy recovered from totatl energy. ##############      
        if self.dispersion is not None:
            van_der_waals_energy = section['x_onetep_total_dispersion_corrected_energy']
            total_energy = section['energy_total']
            for i in range(len(total_energy)):
                self.disp_energy = abs(van_der_waals_energy[i] - total_energy[i])
            backend.addValue('energy_van_der_Waals', self.disp_energy)
        

        finite_basis_corr_energy = section['x_onetep_total_energy_corrected_for_finite_basis_store'] ###Conversion to Jule
        J_converter = 1.602176565e-19
        J_float = float(J_converter)

        if finite_basis_corr_energy:
            finite_basis_corr_energy[0] = float(finite_basis_corr_energy[0]) * J_float
            backend.addValue('x_onetep_total_energy_corrected_for_finite_basis', finite_basis_corr_energy[0])
        
        extFile = ".md"       # Find the file with extension .cell
        dirName = os.path.dirname(os.path.abspath(self.fName))
        cFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                cFile = file
        fName = os.path.normpath(os.path.join(dirName, cFile))

         # parsing *.cell file to get the k path segments
        if file.endswith(extFile):   
            pass
        else:    
            backend.addValue('number_of_scf_iterations', len(self.energy_total_scf_iteration_list))
            # backend.addArrayValues('stress_tensor',np.asarray(self.stress_tensor_value))
        
        if self.local_pseudo:
            backend.addValue('x_onetep_pseudo_local',self.local_pseudo[0])
        if self.kinetic:
            backend.addValue('electronic_kinetic_energy',self.kinetic[0])
        if self.nonlocal_pseudo:
            backend.addValue('x_onetep_pseudo_non_local',self.nonlocal_pseudo[0])
        if self.Hartree:
            backend.addValue('energy_correction_hartree',self.Hartree[0])
        if self.Exchangecorrelation:
            backend.addValue('energy_XC',self.Exchangecorrelation[0])
        if self.Ewald:
            backend.addValue('x_onetep_ewald_correction',self.Ewald[0])
        if self.DispersionCorrection:
            backend.addValue('x_onetep_dispersion_correction',self.DispersionCorrection[0])
        if self.Integrateddensity:
            backend.addValue('x_onetep_integrated_density',self.Integrateddensity[0])
    
    def onClose_section_scf_iteration(self, backend, gIndex, section):
        """trigger called when _section_scf_iteration is closed"""
        # get cached values for energy_total_scf_iteration
        ev = section['energy_total_scf_iteration']
        if ev:
            self.scfIterNr = len(ev)
            self.energy_total_scf_iteration_list.append(ev)
   


    def onClose_x_onetep_section_SCF_iteration_frame(self, backend, gIndex, section):
        
        self.frame_free_energy = section['x_onetep_frame_energy_free']
        self.frame_energies = section['x_onetep_SCF_frame_energy']
        self.frame_energies_gain = section['x_onetep_SCF_frame_energy_gain']               
        self.frame_T0 = section ['x_onetep_frame_energy_total_T0']
        self.wall_time_store = section ['x_onetep_frame_time_scf_iteration_wall_end']

        frame_time = section['x_onetep_frame_time']
        
        
        if self.frame_energies:
 
            for i in range(len(self.frame_energies)):
                
                J_converter = float(1.602176565e-19)

                self.frame_energies[i]=self.frame_energies[i].split()
                self.frame_energies[i]=[float(j) for j in self.frame_energies[i]]
                
                self.frame_energies_gain[i]=self.frame_energies_gain[i].split()
                self.frame_energies_gain[i]=[float(j) for j in self.frame_energies_gain[i]]
                
                self.wall_time_store[i] = self.wall_time_store[i].split()              
                self.wall_time_store[i]=[float(j) for j in self.wall_time_store[i]]
                wall_times = self.wall_time_store[i]

                energies = self.frame_energies[i] ###Conversion to Jule
                energies = [x * J_converter for x in energies]

                energies_gain = self.frame_energies_gain[i] ###Conversion to Jule
                energies_gain = [x * J_converter for x in energies_gain]
                             
                self.energy_frame.extend(energies)   
                
                self.energy_frame_gain.extend(energies_gain) 
                
                self.wall_time_end.extend(wall_times)
            
            free_energies = self.frame_free_energy
            free_energies = [x * J_converter for x in free_energies]
            self.energy_frame_free.extend(free_energies)
           
            t0_energies = self.frame_T0
            t0_energies = [x * J_converter for x in t0_energies]
            self.energy_frame_T0.extend(t0_energies)
            

            if frame_time:
                for i in range(len(frame_time)):                   
                    self.time_0.extend(frame_time)                
                   
        #get cached values of onetep_store_atom_forces
# Recover SCF k points and eigenvalue from *.band file (ONLY FOR SINGLE POINT CALCULATIONS AT THIS STAGE)
    def onClose_x_onetep_section_collect_scf_eigenvalues(self, backend, gIndex, section):

        bandSuperContext = OnetepBandParser.OnetepBandParserContext(False)
        bandParser = AncillaryParser(
            fileDescription = OnetepBandParser.build_onetepBandFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = OnetepBandParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = bandSuperContext)
        
        extFile = ".bands"       # ".band_sp" = spin polarised, ".band" = not spin polarised (ONLY FOR TEST FILES IN /test/examples)
        dirName = os.path.dirname(os.path.abspath(self.fName))
        bFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                bFile = file
                fName = os.path.normpath(os.path.join(dirName, bFile))

                with open(fName) as fIn:    
                    bandParser.parseFile(fIn)  # parsing *.band file to get SCF eigenvalues and relative k points


                self.k_nr_scf     = bandSuperContext.k_nr
                self.n_spin_channels = bandSuperContext.n_spin
            
                self.e_nr_scf     = bandSuperContext.e_nr
                self.k_points_scf = bandSuperContext.eigenvalues_kpoints
                self.e_spin_1     = bandSuperContext.e_spin_1
                self.e_spin_2     = bandSuperContext.e_spin_2
                
            else: 
                pass # if no .bands is found in the same folder as .onetep file than skip and continue     
                

        #self.n_spin = bandSuperContext.n_spin

    def onClose_x_onetep_section_stress_tensor(self, backend, gIndex, section): 
     #get cached values for stress tensor
        stress_tens =[]
        stress_tens = section['x_onetep_store_stress_tensor']
    
        for i in range(len(stress_tens)):
            stress_tens[i] = stress_tens[i].split()
            stress_tens[i] = [float(j) for j in stress_tens[i]]
            stress_tens_int = stress_tens[i]
            stress_tens_int = [x / 10e9 for x in stress_tens_int] #converting GPa in Pa.
            self.stress_tensor_value.append(stress_tens_int)
        self.stress_tensor_value = self.stress_tensor_value[-3:]     

    def onClose_x_onetep_section_raman_tensor(self, backend, gIndex, section): 
     #get cached values for stress tensor
        ram_tens =[]
        ram_tens = section['x_onetep_store_raman_tensor']
        
        for i in range(len(ram_tens)):
            ram_tens[i] = ram_tens[i].split()
            ram_tens[i] = [float(j) for j in ram_tens[i]]
            ram_tens_int = ram_tens[i]
            self.raman_tensor_value.append(ram_tens_int)
        self.raman_tensor_value = self.raman_tensor_value[-3:]         
        if self.raman_tensor_value:
            backend.addArrayValues('x_onetep_ramen_tensor',np.asarray(self.raman_tensor_value))    
    
    

######################################################################################
################ Triggers on closure section_system ######################
######################################################################################

    def onClose_section_system(self, backend, gIndex, section):
    #     number_of_atoms = section['x_onetep_number_of_atoms']
        
       # = number_of_atoms[i].split()
        # self.at_nr = number_of_atoms[i]
     
    #     """trigger called when _section_system is closed"""
    #     cellSuperContext = OnetepCellParser.OnetepCellParserContext(False)
    #     cellParser = AncillaryParser(
    #         fileDescription = OnetepCellParser.build_OnetepCellFileSimpleMatcher(),
    #         parser = self.parser,
    #         cachingLevelForMetaName = OnetepCellParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
    #         superContext = cellSuperContext)

    #     extFile = ".dat"       # Find the file with extension .cell
    #     dirName = os.path.dirname(os.path.abspath(self.fName))
    #     cFile = str()
    #     for file in os.listdir(dirName):
    #         if file.endswith(extFile):
    #             cFile = file
    #         # else:
    #         #     print("ERROR: please provide .dat file to parse info about cell and positions")
    #     fName = os.path.normpath(os.path.join(dirName, cFile))
    #     for file in os.listdir(dirName):
    #         if file.endswith(extFile):
    #             with open(fName) as fIn:
    #                 cellParser.parseFile(fIn)  # parsing *.cell file to get the k path segments
    #             self.cell = cellSuperContext.cell_store  # recover k path segments coordinartes from *.cell file
    #             self.at_nr = cellSuperContext.at_nr
    #             self.atom_labels = cellSuperContext.atom_labels_store
    #             self.onetep_atom_positions = cellSuperContext.onetep_atom_positions_store
    # # Processing the atom positions in fractionary coordinates (as given in the onetep output)
    #             backend.addArrayValues('x_onetep_atom_positions', np.asarray(self.onetep_atom_positions))


    # # Backend add the total number of atoms in the simulation cell
                

    # # Processing the atom labels
    #         #get cached values of onetep_store_atom_labels
                
    #             backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))

    #             backend.addValue('number_of_atoms', self.at_nr)
            
    #             backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:]))
       
       

# Converting the fractional atomic positions (x) to cartesian coordinates (X) ( X = M^-1 x )
        # for i in range(0, self.at_nr):

        #     pos_a = [   self.a[0] * self.onetep_atom_positions[i][0]
        #               + self.b[0] * math.cos(np.deg2rad(self.gamma[0])) * self.onetep_atom_positions[i][1]
        #               + self.c[0] * math.cos(np.deg2rad(self.beta[0])) * self.onetep_atom_positions[i][2],

        #                 self.b[0] * math.sin(self.gamma[0]) * self.onetep_atom_positions[i][1]
        #               + self.c[0] * self.onetep_atom_positions[i][2] * (( math.cos(np.deg2rad(self.alpha[0]))
        #               - math.cos(np.deg2rad(self.beta[0])) * math.cos(np.deg2rad(self.gamma[0])) ) / math.sin(np.deg2rad(self.gamma[0])) ),

        #                (self.volume / (self.a[0]*self.b[0] * math.sin(np.deg2rad(self.gamma[0])))) * self.onetep_atom_positions[i][2] ]

        #     self.atom_positions.append(pos_a)
            
        
        # backend.addValue('x_onetep_cell_volume', self.volume)     
        # backend.addArrayValues('atom_positions', np.asarray(self.atom_positions))
        
#         vel = section['x_onetep_store_atom_ionic_velocities']
#         if vel:
#             for i in range(0, self.at_nr):
#                 vel[i] = vel[i].split()
#                 vel[i] = [float(j) for j in vel[i]]
#                 self.onetep_atom_velocities.append(vel[i])
                
#             backend.addArrayValues('x_onetep_atom_ionic_velocities', np.asarray(self.onetep_atom_velocities))

# # Backend add the simulation cell
        pos_opt = section['x_onetep_store_optimised_atom_positions']
        if pos_opt:

            # backend.addArrayValues('atom_labels', np.asarray(self.atom_labels[-self.at_nr:]))
            
            self.at_nr_opt = len(pos_opt)
            
            for i in range(0, self.at_nr_opt):
                pos_opt[i] = pos_opt[i].split()
                pos_opt[i] = [float(j) for j in pos_opt[i]]
                self.onetep_optimised_atom_positions.append(pos_opt[i])
            backend.addArrayValues('x_onetep_atom_positions', np.asarray(self.onetep_optimised_atom_positions[-self.at_nr_opt:]))
#         # #     print pos_opt[i]    
       
# # Converting the fractional atomic positions (x) to cartesian coordinates (X) ( X = M^-1 x )
#             for i in range(0, self.at_nr_opt):

#                 pos_opt_a = [   self.a[0] * self.onetep_optimised_atom_positions[i][0]
#                       + self.b[0] * math.cos(np.deg2rad(self.gamma[0])) * self.onetep_optimised_atom_positions[i][1]
#                       + self.c[0] * math.cos(np.deg2rad(self.beta[0])) * self.onetep_optimised_atom_positions[i][2],

#                         self.b[0] * math.sin(self.gamma[0]) * self.onetep_optimised_atom_positions[i][1]
#                       + self.c[0] * self.onetep_optimised_atom_positions[i][2] * (( math.cos(np.deg2rad(self.alpha[0]))
#                       - math.cos(np.deg2rad(self.beta[0])) * math.cos(np.deg2rad(self.gamma[0])) ) / math.sin(np.deg2rad(self.gamma[0])) ),

#                        (self.volume / (self.a[0]*self.b[0] * math.sin(np.deg2rad(self.gamma[0])))) * self.onetep_optimised_atom_positions[i][2] ]

#                 self.atom_optim_position.append(pos_opt_a)
            
#             backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:])) 
#             backend.addArrayValues('atom_positions', np.asarray(self.atom_optim_position[-self.at_nr:]))
#             backend.addValue('x_onetep_cell_volume', self.volume) 
#         else:
#             pass
        

    def onClose_x_onetep_section_vibrational_frequencies(self, backend, gIndex, section):
        frequ = section['x_onetep_vibrationl_frequencies_store']
        ##### in this case the name x_onetep_ir_store sands for irreducible representation in the point group.
        irr_rep = section ['x_onetep_ir_store']
        ##### in this case the name x_onetep_ir_intesities_store sands for Infra-red intensities. 
        ir_intensities = section['x_onetep_ir_intensity_store']
     
        self.nr_iter =section['x_onetep_n_iterations_phonons']

        raman_activity = section['x_onetep_raman_activity_store']
       
        #raman_act = section['onetep_raman_active_store']
        
        if frequ and ir_intensities and raman_activity:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.extend(frequ_list)    
                    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.extend(irr_rep_list)
                    
                    ir_intensities[i] = ir_intensities[i].split()
                    ir_intensities[i] = [float(j) for j in ir_intensities[i]]
                    ir_intens_list = ir_intensities[i]
                    self.ir_intens.extend(ir_intens_list) 
                    
                    raman_activity[i] = raman_activity[i].split()
                    raman_activity[i] = [float(j) for j in raman_activity[i]]
                    raman_list = raman_activity[i]
                    self.raman_act.extend(raman_list)  
                
                
                backend.addArrayValues('x_onetep_ir_intensity', np.asarray(self.ir_intens[-len(self.nr_iter):]))
                backend.addArrayValues('x_onetep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):])) 
                backend.addArrayValues('x_onetep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
                backend.addArrayValues('x_onetep_raman_activity', np.asarray(self.raman_act[-len(self.nr_iter):]))
        
        elif frequ and ir_intensities:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.extend(frequ_list)    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.extend(irr_rep_list)
                    
                    ir_intensities[i] = ir_intensities[i].split()
                    ir_intensities[i] = [float(j) for j in ir_intensities[i]]
                    ir_intens_list = ir_intensities[i]
                    self.ir_intens.extend(ir_intens_list)
                  
                backend.addArrayValues('x_onetep_ir_intensity', np.asarray(self.ir_intens[-len(self.nr_iter):]))
                backend.addArrayValues('x_onetep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):]))    
                backend.addArrayValues('x_onetep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
        
        elif frequ:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.extend(frequ_list)    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.extend(irr_rep_list) 
                
                backend.addArrayValues('x_onetep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):])) 
                backend.addArrayValues('x_onetep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
            
    def onClose_x_onetep_section_population_analysis(self, backend, gIndex, section):
        orb_contr = section['x_onetep_total_orbital_store']
        tot_charge = section ['x_onetep_mulliken_charge_store']
        spin = section['x_onetep_spin_store']
        if orb_contr:
            self.total_contribution.extend(orb_contr)
   
            backend.addArrayValues('x_onetep_total_orbital',np.asarray(self.total_contribution))
        
        if tot_charge:
            self.total_charge.extend(tot_charge)
            backend.addArrayValues('x_onetep_mulliken_charge', np.asarray(self.total_charge))

        if spin:
            self.total_spin_mulliken.extend(spin)
   
            backend.addArrayValues('x_onetep_spin',np.asarray(self.total_spin_mulliken))
######################################################################################
###################### Storing k band points and band energies #######################
############################# FIRST SPIN CHANNEL #####################################
######################################################################################

# Storing the k point coordinates (SPIN 1)
    def onClose_x_onetep_section_k_points(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
# Processing k points (given in fractional coordinates)
        #get cached values of onetep_store_k_points
        k_st = section['x_onetep_store_k_points']
        self.k_count = len(k_st)
        self.k_nr   += 1
        for i in range(0, self.k_count):
            k_st[i] = k_st[i].split()
            k_st[i] = [float(j) for j in k_st[i]]
            k_st_int = k_st[i]
            self.onetep_band_kpoints.append(k_st_int)
          
        self.n_spin_channels_bands = 1    
        
# Storing the eigenvalues (SPIN 1)
    def onClose_x_onetep_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        #get cached values of onetep_store_k_points
        e_st = section['x_onetep_store_eigenvalues']
        self.e_nr = len(e_st)
        self.onetep_band_energies.append(e_st)


######################################################################################
###################### Storing k band points and band energies #######################
############################# SECOND SPIN CHANNEL ####################################
######################################################################################

# Storing the k point coordinates (SPIN 2)
    def onClose_x_onetep_section_k_points_1(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
# Processing k points (given in fractional coordinates)
        #get cached values of onetep_store_k_points
        k_st_1 = section['x_onetep_store_k_points_1']
        self.k_count_1 = len(k_st_1)
        self.k_nr_1   += 1
        for i in range(0, self.k_count_1):
            k_st_1[i] = k_st_1[i].split()
            k_st_1[i] = [float(j) for j in k_st_1[i]]
            k_st_1_int = k_st_1[i]
            self.onetep_band_kpoints_1.append(k_st_1_int)

        self.k_nr_1 = self.k_nr  # clean double counting
        
        self.n_spin_channels_bands = 2

# Storing the eigenvalues (SPIN 2)
    def onClose_x_onetep_section_eigenvalues_1(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        #get cached values of onetep_store_k_points
        e_st_1 = section['x_onetep_store_eigenvalues_1']
        self.e_nr_1 = len(e_st_1)
        self.onetep_band_energies_1.append(e_st_1)

        self.e_nr_1 = self.e_nr
        

######################################################################################
########### BAND STRUCTURE ###########################################################
######################################################################################

    # def onClose_section_k_band(self, backend, gIndex, section):
    #     """Trigger called when section_k_band is closed.

    #        The k path coordinates are extracted from the *.cell input file.
    #     """

    #     cellSuperContext = OnetepCellParser.onetepCellParserContext(False)
    #     cellParser = AncillaryParser(
    #         fileDescription = OnetepCellParser.build_onetepCellFileSimpleMatcher(),
    #         parser = self.parser,
    #         cachingLevelForMetaName = OnetepCellParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
    #         superContext = cellSuperContext)

    #     extFile = ".cell"       # Find the file with extension .cell
    #     dirName = os.path.dirname(os.path.abspath(self.fName))
    #     cFile = str()
    #     for file in os.listdir(dirName):
    #         if file.endswith(extFile):
    #             cFile = file
    #     fName = os.path.normpath(os.path.join(dirName, cFile))

    #     with open(fName) as fIn:
    #         cellParser.parseFile(fIn)  # parsing *.cell file to get the k path segments
           
    #     self.k_start_end = cellSuperContext.k_sgt_start_end  # recover k path segments coordinartes from *.cell file
    #     self.k_path_nr = len(self.k_start_end)
    #     # backend.openSection('section_single_configuration_to_calculation_ref')
    #     if self.n_spin_channels_bands:
    #         backend.openSection('section_method')
    #         backend.addValue('number_of_spin_channels',self.n_spin_channels_bands)     
    #         backend.closeSection('section_method',gIndex+1)
        
        ########################################################################################
        def get_last_index(el, check):  # function that returns end index for each k path segment
            found = None
            for i, next in enumerate(check):
                if next == el:
                    found = i

            assert found != None
            return found
        ########################################################################################
        if self.k_start_end: 
            if self.onetep_band_energies_1 != []:  # handling k band energies
                for i in range(self.k_nr):
                    a = [ self.onetep_band_energies[i], self.onetep_band_energies_1[i] ]  # spin polarised
                    self.band_en.append(a)
            else:
                self.band_en = self.onetep_band_energies  # single spin

    
            path_end_index = []
            for i in range(self.k_path_nr):
                boundary = self.k_start_end[i][1]
                a = get_last_index(boundary, self.onetep_band_kpoints)
                path_end_index.append(a)

            path_end_index = [0] + path_end_index  # list storing the end index of each k segment


            k_point_path = []
            for i in range(self.k_path_nr):
                a = self.onetep_band_kpoints[ path_end_index[i] : path_end_index[i+1]+1 ]
                k_point_path.append(a)          # storing the k point fractional coordinates for each segment


            band_en_path = []
            for i in range(self.k_path_nr):
                a = self.band_en[ path_end_index[i] : path_end_index[i+1]+1 ]
                band_en_path.append(a)          # storing the band energies for each segment, k point and spin channel
           
            for i in range(self.k_path_nr):    
              
                backend.openSection('section_k_band_segment')
                backend.addArrayValues('band_k_points', np.asarray(k_point_path[i]))
                backend.addArrayValues('band_energies', np.asarray(band_en_path[i]))
                backend.addArrayValues('band_segm_start_end', np.asarray(self.k_start_end[i])) 
                backend.addValue('number_of_k_points_per_segment', len(k_point_path[i])) 
                backend.closeSection('section_k_band_segment',i)
            
        else: 
            pass    

    def onClose_x_onetep_section_energy_components(self, backend, gIndex, section):        
        ######### Caching energy components #######
        self.kinetic = section['x_onetep_electronic_kinetic_energy']
        self.local_pseudo = section['x_onetep_pseudo_local_store'] 
        self.nonlocal_pseudo = section['x_onetep_pseudo_non_local_store']
        
        self.Hartree = section['x_onetep_energy_correction_hartree_store']
        self.Exchangecorrelation = section['x_onetep_energy_XC_store'] 
        self.Ewald = section['x_onetep_ewald_correction_store']
        self.DispersionCorrection =section['x_onetep_dispersion_correction_store'] 
        self.Integrateddensity = section['x_onetep_integrated_density_store']

    
    def onClose_section_run(self, backend, gIndex, section):
        
        f_st_band = section['x_onetep_store_atom_forces_band']
       
        # if f_st_band:
        #     gindex_band = 1
        #     for i in range(0, self.at_nr):
        #         f_st_band[i] = f_st_band[i].split()
        #         f_st_band[i] = [float(j) for j in f_st_band[i]]

        #         f_st_int_band = f_st_band[i]
                 
        #         self.atom_forces_band.append(f_st_int_band)
        #         self.atom_forces_band = self.atom_forces_band[-self.at_nr:] 
        #     backend.addArrayValues('x_onetep_atom_forces', np.asarray(self.atom_forces_band))        

        cellSuperContext = OnetepCellParser.OnetepCellParserContext(False)
        cellParser = AncillaryParser(
            fileDescription = OnetepCellParser.build_OnetepCellFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = OnetepCellParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = cellSuperContext)

        extFile = ".dat"       # Find the file with extension .cell
        dirName = os.path.dirname(os.path.abspath(self.fName))
        cFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                cFile = file
            # else:
            #     print("ERROR: please provide .dat file to parse info about cell and positions")
        fName = os.path.normpath(os.path.join(dirName, cFile))
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                with open(fName) as fIn:
                    cellParser.parseFile(fIn)  # parsing *.cell file to get the k path segments
                self.cell = cellSuperContext.cell_store  # recover k path segments coordinartes from *.cell file
                self.at_nr = cellSuperContext.at_nr
                self.atom_labels = cellSuperContext.atom_labels_store
                self.onetep_atom_positions = cellSuperContext.onetep_atom_positions_store
    # Processing the atom positions in fractionary coordinates (as given in the onetep output)
                i = backend.openSection('section_system')
               
                backend.addArrayValues('x_onetep_atom_positions', np.asarray(self.onetep_atom_positions))

                backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))

                backend.addValue('number_of_atoms', len(self.atom_labels))
            
                backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:]))    
               
                backend.closeSection('section_system', i )
        
        
        # time_list = self.time_0
        if section['x_onetep_geom_converged'] is not None:
            if section['x_onetep_geom_converged'][-1] == 'successfully':
                self.geoConvergence = True
            else:
                self.geoConvergence = False
        if self.geoConvergence is not None:
            backend.openSection('section_frame_sequence')
            backend.addValue('geometry_optimization_converged', self.geoConvergence)        
            backend.closeSection('section_frame_sequence',gIndex)
        

        MDSuperContext = OnetepMDParser.OnetepMDParserContext(False)
        MDParser = AncillaryParser(
            fileDescription = OnetepMDParser.build_OnetepMDFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = OnetepMDParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = MDSuperContext)

        extFile = ".md"       # Find the file with extension .cell
        dirName = os.path.dirname(os.path.abspath(self.fName))
        cFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                cFile = file
        fName = os.path.normpath(os.path.join(dirName, cFile))

         # parsing *.cell file to get the k path segments
        if file.endswith(extFile):   
            with open(fName) as fIn:
                MDParser.parseFile(fIn)
            
            self.frame_temp = MDSuperContext.frame_temperature  # recover k path segments coordinartes from *.cell file
            self.frame_press= MDSuperContext.frame_pressure
            self.frame_kinetic = MDSuperContext.kinetic
            self.frame_potential = MDSuperContext.hamiltonian
            self.frame_atom_forces = MDSuperContext.total_forces
            self.frame_atom_veloc = MDSuperContext.total_velocities
            self.frame_stress_tensor = MDSuperContext.frame_stress_tensor
            self.frame_position =MDSuperContext.total_positions
            self.frame_cell=MDSuperContext.frame_cell
            self.frame_vet_velocities =MDSuperContext.vector_velocities
            

            
        #     if self.frame_temp:
                
        #         gIndexGroupscf = backend.openSection('section_scf_iteration')
        #         for i in range(len(self.frame_atom_forces)):
                    
        #             backend.openSection('section_system')
        #             backend.addArrayValues('atom_velocities', np.asarray(self.frame_atom_veloc[i]))
        #             backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))
        #             backend.addArrayValues('atom_positions', np.asarray(self.frame_position[i]))
        #             backend.addArrayValues('simulation_cell', np.asarray(self.frame_cell[i]))
        #             backend.addArrayValues('x_onetep_velocities_cell_vector',np.asarray(self.frame_vet_velocities[i]))
        #             backend.closeSection('section_system',i+2)
                    
        #             backend.openSection('section_single_configuration_calculation')
        #             backend.addArrayValues('atom_forces', np.asarray(self.frame_atom_forces[i]))
        #             backend.addArrayValues('stress_tensor',np.asarray(self.frame_stress_tensor[i]))
        #             backend.addValue('number_of_scf_iterations', len(self.frame_energies))
        #             if i > 0:
                     
        #                 backend.addValue('energy_free', self.energy_frame_free[i-1]) 
        #                 backend.addValue('energy_total_T0',self.energy_frame_T0[i-1])
                    
        #             if i > 0:    
        #                 for j in range(len(self.frame_energies)):
        #                     s = j + i*len(self.frame_energies) - len(self.frame_energies)
                            
        #                     backend.openSection('section_scf_iteration')
        #                     backend.addValue('energy_total_scf_iteration', self.energy_frame[s])
        #                     backend.addValue('energy_change_scf_iteration', self.energy_frame_gain[s])
        #                     backend.addValue('time_scf_iteration_wall_end',  self.wall_time_end[s])
        #                     backend.closeSection('section_scf_iteration',s+gIndexGroupscf)
                                                         
        #             backend.closeSection('section_single_configuration_calculation',i+1) 

        #         backend.openSection('section_frame_sequence')
        #         backend.addValue('number_of_frames_in_sequence',(len(self.frame_potential)))
        #         backend.addArrayValues('frame_sequence_temperature', np.asarray(self.frame_temp))
        #         backend.addArrayValues('frame_sequence_pressure', np.asarray(self.frame_press))
        #         backend.addArrayValues('frame_sequence_kinetic_energy', np.asarray(self.frame_kinetic))
        #         backend.addArrayValues('frame_sequence_potential_energy', np.asarray(self.frame_potential))
        #         backend.addArrayValues('frame_sequence_time', np.asarray(time_list))
        #         backend.closeSection('section_frame_sequence',gIndex)
            
        #     else:
        #         pass
        
        # else:
        #     pass            

        
        # TSSuperContext = OnetepTSParser.OnetepTSParserContext(False)
        # TSParser = AncillaryParser(
        #     fileDescription = OnetepTSParser.build_onetepTSFileSimpleMatcher(),
        #     parser = self.parser,
        #     cachingLevelForMetaName = OnetepTSParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
        #     superContext = TSSuperContext)

        # extFile = ".ts"       # Find the file with extension .cell
        # dirName = os.path.dirname(os.path.abspath(self.fName))
        # cFile = str()
        # for file in os.listdir(dirName):
        #     if file.endswith(extFile):
        #         cFile = file
        #     fName = os.path.normpath(os.path.join(dirName, cFile))
        #     if file.endswith(".ts"):
        #         with open(fName) as fIn:
        #             TSParser.parseFile(fIn)          

        #             self.ts_total_energy = TSSuperContext.total_energy    
        #             self.ts_cell_vector = TSSuperContext.frame_cell
        #             self.ts_forces = TSSuperContext.total_forces
        #             self.ts_position = TSSuperContext.total_positions
                   
        #             self.ts_path = TSSuperContext.path_ts
                  
        #             self.ts_total_energy_f = TSSuperContext.total_energy_final
        #             self.ts_forces_f = TSSuperContext.md_forces_final
        #             self.ts_cell_vector_f = TSSuperContext.cell_final
        #             self.ts_positions_f = TSSuperContext.atomf_position
        #             self.ts_path_f =TSSuperContext.path_final
        #             self.ts_total_energy_p = TSSuperContext.total_energy_pro
        #             self.ts_forces_p = TSSuperContext.md_forces_pro
        #             self.ts_cell_vector_p = TSSuperContext.cell_pro
        #             self.ts_positions_p = TSSuperContext.atomp_position
        #             self.ts_path_p =TSSuperContext.path_pro
        #             for i in range(len(self.ts_total_energy)):
        #                 backend.openSection('x_onetep_section_ts')
        #                 backend.addValue('x_onetep_ts_path', self.ts_path[i]) 
        #                 backend.addValue('x_onetep_ts_energy_total', self.ts_total_energy[i])
        #                 backend.addArrayValues('x_onetep_ts_cell_vectors', np.asarray(self.ts_cell_vector[i]))
        #                 backend.addArrayValues('x_onetep_ts_forces', np.asarray(self.ts_forces[i]))
        #                 backend.addArrayValues('x_onetep_ts_positions', np.asarray(self.ts_position[i]))
                                               
        #                 backend.closeSection('x_onetep_section_ts',i)        

        #             backend.openSection('x_onetep_section_ts_final')
        #             backend.addValue('x_onetep_ts_energy_final', self.ts_total_energy_f)
        #             backend.addArrayValues('x_onetep_ts_cell_vectors_final', np.asarray(self.ts_cell_vector_f))
        #             backend.addArrayValues('x_onetep_ts_positions_final', np.asarray(self.ts_positions_f))
        #             backend.addArrayValues('x_onetep_ts_forces_final', np.asarray(self.ts_forces_f))
                  
        #             backend.addValue('x_onetep_ts_path_ts_final', self.ts_path_f)    
        #             backend.closeSection('x_onetep_section_ts_final',gIndex)    
        #             backend.openSection('x_onetep_section_ts_product')
        #             backend.addValue('x_onetep_ts_energy_product', self.ts_total_energy_p)
        #             backend.addArrayValues('x_onetep_ts_cell_vectors_product', np.asarray(self.ts_cell_vector_p))
        #             backend.addArrayValues('x_onetep_ts_positions_product', np.asarray(self.ts_positions_p))
        #             backend.addArrayValues('x_onetep_ts_forces_product', np.asarray(self.ts_forces_p))
        #             backend.addValue('x_onetep_ts_path_product', self.ts_path_p)    
        #             backend.closeSection('x_onetep_section_ts_product',gIndex)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  MAIN PARSER STARTS HERE  ###############################################################################################################
################################################################################################################################################################
############################################################################ onetep.Parser Version 1.0 #########################################################
################################################################################################################################################################


def build_onetepMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main file of onetep.

    First, several subMatchers are defined, which are then used to piece together
    the final SimpleMatcher.

    Returns:
       SimpleMatcher that parses main file of onetep.
    """

    ########################################
  
    ########################################
    # submatcher for
    scfEigenvaluesSubMatcher = SM(name = 'scfEigenvalues',
       startReStr = r"\stype\sof\scalculation\s*\:\ssingle\spoint\senergy\s*",
       sections = ["x_onetep_section_collect_scf_eigenvalues"]

       ) # CLOSING SM scfEigenvalues


    ########################################
    # submatcher for section method
    


    ########################################
    # submatcher for section_basis_set_cell_dependent

    basisSetCellAssociatedSubMatcher = SM(name = 'planeWaveBasisSet',
        startReStr = r"cutoff_energy+" or r"kernel_cutoff+",
        subFlags = SM.SubFlags.Unordered,
        forwardMatch = True,
        sections = ["section_basis_set_cell_dependent"],
        subMatchers = [
            SM(r"kernel_cutoff\s+(?P<x_onetep_kernel_cutoff>[0-9.]+)"),
            SM(r"cutoff_energy\s+(?P<x_onetep_basis_set_planewave_cutoff>[0-9.]+)"),
            SM(r"kernel_cutoff\s*\:\s*(?P<x_onetep_kernel_cutoff>[0-9.]+)"),
            SM(r"cutoff_energy\s*\:\s*(?P<x_onetep_basis_set_planewave_cutoff>[0-9.]+)"),
            SM(r"kernel_cutoff\:\s*(?P<x_onetep_kernel_cutoff>[0-9.]+)"),
            SM(r"cutoff_energy\:\s*(?P<x_onetep_basis_set_planewave_cutoff>[0-9.]+)"),
           
                      ]) # CLOSING SM planeWaveBasisSet

    calculationMethodSubMatcher = SM(name = 'calculationMethods',
        startReStr = r"xc\_functional+",
        forwardMatch = True,
        sections = ["section_method"],
        subMatchers = [

            SM(name = "onetepXC",
               startReStr = r"xc\_functional+",
               subFlags = SM.SubFlags.Unordered,
               forwardMatch = True,
               sections = ["x_onetep_section_functionals"],
               subMatchers = [
                 SM(r"xc\_functional\s*\:\s(?P<x_onetep_functional_name>[A-Za-z0-9()]+)"),
                 SM(r"xc\_functional\:\s*(?P<x_onetep_functional_name>[A-Za-z0-9()]+)"),
                 SM(r"xc\_functional\s*(?P<x_onetep_functional_name>[A-Za-z0-9()]+)") ,
                 SM(r"xc\_functional\s*\:\s*(?P<x_onetep_functional_name>[A-Za-z0-9()]+)"),
                 
                             ]), # CLOSING onetep_section_functionals
            ])
    #        # SM(name = "onetepXC_definition",
    #        #    startReStr = r"\susing custom XC functional definition\:",
    #        #    #endReStr = r"\srelativistic treatment\s*\:\s*",
    #        #    forwardMatch = True,
    #        #    sections = ["x_onetep_section_functional_definition"],
    #        #    subMatchers = [     
    #        #       SM(r"\s*(?P<x_onetep_functional_type>[A-Za-z0-9]+)\s*(?P<x_onetep_functional_weight>[0-9.]+)",
    #        #              repeats = True),
    #        #       #SM(r"\srelativistic treatment\s*\: *(?P<onetep_relativity_treatment_scf> [A-Za-z0-9() -]*)")              
    #        #                    ]), # CLOSING onetep_section_functional_definition

    #        # SM(name = "onetep_relativ",
    #        #    startReStr = r"\srelativistic treatment\s*\:\s*",
    #        #    forwardMatch = True,
    #        #    sections = ["x_onetep_section_relativity_treatment"],
    #        #    subMatchers = [
    #        #       SM(r"\srelativistic treatment\s*\: *(?P<x_onetep_relativity_treatment_scf> [A-Za-z0-9() -]+)")
    #        #                   ]), # CLOSING onetep_section_relativistic_treatment           

    #        #  SM(name = "van der Waals",
    #        #     startReStr = r"\sDFT\+D: Semi-empirical dispersion correction\s*\:\s*",
    #        #     #forwardMatch = True,
    #        #     #sections = ["onetep_section_relativity_treatment"],
    #        #     subMatchers = [
    #        #       SM(r"\sSEDC with\s*\: *(?P<van_der_Waals_method> [A-Za-z0-9() -]+)"),
    #        #                   ]), # CLOSING van der waals
            
            
              

    #           ]), # CLOSING SM calculationMethods

    phononCalculationSubMatcher = SM(name = 'phonon_calculation',
        sections = ["x_onetep_section_phonons"],
        startReStr = r"\s\*\*\** Phonon Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\sphonon calculation method\s*\:\s*(?P<x_onetep_phonon_method>[a-zA-Z]+\s[a-zA-Z]+)"),
            SM(r"\sphonon convergence tolerance\s*\:\s*(?P<x_onetep_phonon_tolerance>[-+0-9.eEd]+)"),
            SM(r"\smax\. number of phonon cycles\s*\:\s*(?P<x_onetep_phonon_cycles>[0-9.]+)"),
            SM(r"\sDFPT solver method\s*\:\s*(?P<x_onetep_DFPT_solver_method>[a-zA-Z0-9.() ]+)"),
            SM(r"\sband convergence tolerance\s*\:\s*(?P<x_onetep_band_tolerance>[-+0-9.eEd]+)"),
                          ])

    GeomOptimParameterSubMatcher = SM(name = 'optimistation_parameters',
        sections = ["section_sampling_method"],
        # sections = ["x_onetep_section_geom_optimisation_method"],
        startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\ Geometry Optimization Parameters \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
        subMatchers = [
            
            SM(r"\soptimization method\s*\:\s*(?P<geometry_optimization_method>[a-zA-Z ()]+)"),
            SM(r"\stotal energy convergence tolerance\s*\:\s*(?P<geometry_optimization_energy_change>[-+0-9.eEd]+)"),
            SM(r"\smax. number of steps\s*\:\s*(?P<x_onetep_max_number_of_steps>[0-9.]+)"),
            SM(r"\s*max ionic \|force\| tolerance\s*\:\s*(?P<geometry_optimization_threshold_force>[-+0-9.eEd]+)"),
            SM(r"\smax ionic \|displacement\| tolerance\s*\:\s*(?P<geometry_optimization_geometry_change>[-+0-9.eEd]+)"),
            SM(r"\s*max \|stress component\| tolerance\s*\:\s*(?P<x_onetep_geometry_stress_com_tolerance>[-+0-9.eEd]+)"),
                                            ]) # CLOSING converged optmisation
    
    ElectronicParameterSubMatcher = SM(name = 'Elec_parameters' ,            
        sections = ["section_system"],
        startReStr = r"\-*\sAtom counting information\s\-*\s*",
        subMatchers = [
            
            SM(r"Totals\:\s*(?P<x_onetep_number_of_atoms>[0-9]+)\s*(?P<x_onetep_number_of_ngwf>[0-9]+)\s*(?P<x_onetep_number_of_projectors>[0-9]+)"),
            # SM(r"\s*net charge of system\s*\:\s*(?P<x_onetep_net_charge>[+0-9.eEd]+)"),
            # SM(r"\s*number of bands\s*\:\s*(?P<x_onetep_number_of_bands>[0-9.]+)"),
            ])
    
    ElectronicMinimisParameterSubMatcher = SM(name = 'Elec_min_parameters' ,            
        sections = ["x_onetep_section_scf_parameters"],
        startReStr = r"\s\*\*\** Electronic Minimization Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*total energy \/ atom convergence tol.\s*\:\s*(?P<x_onetep_energy_threshold_store>[-+0-9.eEd]+)"),
            SM(r"\s*max. number of SCF cycles\s*\:\s*(?P<x_onetep_max_iter_store>[0-9.]+)"),
            SM(r"\s*smearing scheme\s*\:\s*(?P<x_onetep_smearing_kind>[A-Za-z]+)"),
            SM(r"\s*smearing width\s*\:\s*(?P<x_onetep_smearing_width>[0-9.]+)"),
            ])
    
    DensityMixingParameterSubMatcher = SM(name = 'Density_mixing' ,            
        sections = ["x_onetep_section_density_mixing_parameters"],
        startReStr = r"\s\*\*\** Density Mixing Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*density-mixing scheme\s*\:\s*(?P<x_onetep_density_mixing_scheme>[A-Za-z]+)"),
            SM(r"\s*max\. length of mixing history\s*\:\s*(?P<x_onetep_density_mixing_length>[0-9.]+)"),
            SM(r"\s*charge density mixing amplitude\s*\:\s*(?P<x_onetep_charge_density_mixing_amplitude>[0-9.]+)"),
            SM(r"\s*cut\-off energy for mixing\s*\:\s*(?P<x_onetep_cut_off_energy_for_mixing>[0-9.]+)"),
            ])
    PopulationAnalysisParameterSubMatcher = SM(name = 'Pop_analysis' ,            
        sections = ["x_onetep_section_population_analysis_parameters"],
        startReStr = r"\s\*\*\** Population Analysis Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*Population analysis with cutoff\s*\:\s*(?P<x_onetep_population_analysis_cutoff>[-+0-9.eEd]+)"),
           
            ])
    MDParameterSubMatcher = SM(name = 'MD_parameters' ,            
        sections = ["section_sampling_method"],
        startReStr = r"\s\*\*\** Molecular Dynamics Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*ensemble\s*\:\s*(?P<ensemble_type>[A-Za-z]+)"),
            SM(r"\s*temperature\s*\:\s*(?P<x_onetep_thermostat_target_temperature>[0-9.]+\.)"),
            SM(r"\s*pressure\s*\:\s*(?P<x_onetep_frame_pressure>[A-Za-z0-9]+\.)"),
            SM(r"\s*using\s*\:\s*(?P<x_onetep_barostat_type>[A-Za-z]+\-[A-Za-z]+)\s+"),
            SM(r"\s*with characteristic cell time\s*\:\s*(?P<x_onetep_barostat_tau>[-+0-9.eEd]+)"),
            SM(r"\s*using\s*\:\s*(?P<x_onetep_thermostat_type>[A-Za-z]+\-[A-Za-z]+)\s+"),
            SM(r"\s*with characteristic ionic time\s*\:\s*(?P<x_onetep_thermostat_tau>[-+0-9.eEd]+)"),
            SM(r"\s*time\sstep\s*\:\s*(?P<x_onetep_integrator_dt>[-+0-9.eEdD]+)"),
            SM(r"\s*number of MD steps\s*\:\s*(?P<x_onetep_number_of_steps_requested>[0-9]+)"),
            SM(r"\s*MD SCF energy \/ atom convergence tol\.\s*\:\s*(?P<x_onetep_frame_energy_tolerance>[-+0-9.eEd]+)"),
            SM(r"\s*MD SCF eigenenergies tolerance\s*\:\s*(?P<x_onetep_frame_eigen_tolerance>[-+0-9.eEd]+)"),

            ])
    TSParameterSubMatcher = SM(name = 'TS_parameters' ,            
        sections = ["x_onetep_section_ts_parameters"],
        startReStr = r"\s\*\*\** Transition State Search Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*search method\s*\:\s*(?P<x_onetep_ts_method>[A-Za-z]+)"),
            SM(r"\s*[A-Z]+ protocol\s*\:\s*(?P<x_onetep_ts_protocol>[A-Za-z]+)"),
 
            SM(r"\s*max\. number of QST iterations\s*\:\s*(?P<x_onetep_ts_number_qst>[-+0-9.eEd]+)"),
            SM(r"\s*max\. number of CG iterations\s*\:\s*(?P<x_onetep_ts_number_cg>[-+0-9.eEd]+)"),
            SM(r"\s*max\. ionic \|force\| tolerance\s*\:\s*(?P<x_onetep_ts_force_tolerance>[-+0-9.eEd]+)"),
            SM(r"\s*max\. ionic \|displacement\| tolerance\s*\:\s*(?P<x_onetep_ts_displacement_tolerance>[-+0-9.eEd]+)"),
            ])
    OpticsParameterSubMatcher = SM(name = 'optics_parameters' ,            
        sections = ["x_onetep_section_optics_parameters"],
        startReStr = r"\s\*\*\** Optics Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*search method\s*\:\s*(?P<x_onetep_optics_n_bands>[0-9.]+)"),
            SM(r"\s*band convergence tolerance\s*\:\s*(?P<x_onetep_optics_tolerance>[-+0-9.eEd]+)"),
 
            
            ])
    ElecSpecParameterSubMatcher = SM(name = 'electronic_spectroscopy_parameters' ,            
        sections = ["x_onetep_section_electronic_spectroscpy_parameters"],
        startReStr = r"\s\*\*\** Electronic Spectroscopy Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*electronic spectroscopy with theory level\s*\:\s*(?P<x_onetep_theory_level>[A-Z]+)"),
            SM(r"\s*spectroscopy calculation\s*\:\s*(?P<x_onetep_spectroscopy_calculation>[A-Za-z0-9 + A-Za-z0-9 ]+)"),
            SM(r"\s*max\. number of iterations \s*\:\s*(?P<x_onetep_spec_max_iter>[0-9.]+)"),
            SM(r"\s*max\. steps per iteration\s*\:\s*(?P<x_onetep_spec_max_steps>[0-9.]+)"),
            SM(r"\s*number of bands \/ k-point\s*\:\s*(?P<x_onetep_spec_max_bands>[0-9.]+)"),
            SM(r"\s*band convergence tolerance\s*\:\s*(?P<x_onetep_spec_tolerance>[-+0-9.eEd]+)"),
            ])
    

    TDDFTParameterSubMatcher = SM(name = 'TDDFT' ,            
        sections = ["x_onetep_section_tddft_parameters"],
        startReStr = r"\s\*\*\** Time-Dependent DFT Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*number of excited states \s*\:\s*(?P<x_onetep_tddft_n_excited_states>[0-9.]+)"),
            SM(r"\s*state selected for calculation of forces\s*\:\s*(?P<x_onetep_tddft_n_states_forces>[0-9.]+)"),
            SM(r"\s*state convergence tolerance\s*\:\s*(?P<x_onetep_tddft_state_tolerance>[-+0-9.eEd]+)"),
            SM(r"\s*convergence tolerance window\s*\:\s*(?P<x_onetep_tddft_state_tolerance_window>[0-9.]+)"),
            SM(r"\s*max. number of iterations\s*\:\s*(?P<x_onetep_tddft_max_iter>[0-9.]+)"),
            SM(r"\s*no. of extra (convergence indifferent) states\s*\:\s*(?P<x_onetep_tddft_extra_states>[0-9.]+)"),
            SM(r"\s*using tddft functional\s*\:\s*(?P<x_onetep_tddft_functional>[A-Za-z0-9 + A-Za-z0-9 ]+)"),
            SM(r"\s*Time-Dependent DFT method\s*\:\s*(?P<x_onetep_tddft_method>[A-Za-z]+)"),
            
            SM(r"\s*matrix eigenvalue method\s*\:\s*(?P<x_onetep_tddft_eigenmethod>[A-Za-z]+)"),
            SM(r"\s*Time-Dependent DFT approximation\s*\:\s*(?P<x_onetep_tddft_approximation>[A-Za-z]+)"),
            SM(r"\s*Time-Dependent DFT position operator\s*\:\s*(?P<x_onetep_tddft_position_op>[A-Za-z]+)"),
            ])

    bandStructureSubMatcher = SM (name = 'BandStructure',
        startReStr = r"\s*\+\s*B A N D   S T R U C T U R E   C A L C U L A T I O N\s*",
        #startReStr = r"\stype\sof\scalculation\s*\:\sband\sstructure\s*",
        # sections = ['section_k_band'],
        sections = ['section_single_configuration_calculation'],
        subMatchers = [

         SM(startReStr = "\s*\+\s*Spin\=1\skpt\=\s*",
               sections = ["section_k_band"],
               # forwardMatch = True,
               subMatchers = [

            # First spin channel
            SM(startReStr = "\s*\+\s*Spin\=1\skpt\=\s*",
               sections = ["x_onetep_section_k_band"],
               forwardMatch = True,
               subMatchers = [

                  SM(startReStr = "\s*\+\s*Spin\=1\skpt\=\s*",
                     sections = ["x_onetep_section_k_points"],
                     forwardMatch = True,
                     repeats = True,
                     subMatchers = [
                        # Matching k points
                        SM(r"\s*\+\s*Spin\=1\s*kpt\=\s*[0-9]+\s*\((?P<x_onetep_store_k_points>\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\)\s*",
                           repeats = True),

                           SM(name = 'Eigen_1',
                              startReStr = r"\s*\+\s*\+\s*",
                              sections = ['x_onetep_section_eigenvalues'],
                              repeats = True,
                              subMatchers = [
                                 # Matching eigenvalues
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<x_onetep_store_eigenvalues>\s+[-+0-9.eEd]+)",
                                    repeats = True)

                                             ]), # CLOSING x_onetep_section_eigenvalues


                                    ]), # CLOSING onetep_section_k_points


                              ]), # CLOSING 1st section_eigenvalues

            # Second spin channel
            SM(startReStr = "\s*\+\s*Spin\=2\skpt\=\s*",
               sections = ["x_onetep_section_k_band"],
               forwardMatch = True,
               subMatchers = [

                  SM(startReStr = "\s*\+\s*Spin\=2\skpt\=\s*",
                     sections = ["x_onetep_section_k_points_1"],
                     forwardMatch = True,
                     repeats = True,
                     subMatchers = [
                        # Matching k points
                        SM(r"\s*\+\s*Spin\=2\s*kpt\=\s*[0-9]+\s*\((?P<x_onetep_store_k_points_1>\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\)\s*",
                           repeats = True),

                           SM(name = 'Eigen_2',
                              startReStr = r"\s*\+\s*\+\s*",
                              sections = ['x_onetep_section_eigenvalues_1'],
                              repeats = True,
                              subMatchers = [
                                 # Matching eigenvalues
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<x_onetep_store_eigenvalues_1>\s+[-+0-9.eEd]+)",
                                    repeats = True)

                                             ]), # CLOSING onetep_section_eigenvalues_1


                                    ]), # CLOSING onetep_section_k_points_1


                              ]), # CLOSING 2nd section_eigenvalues


        ]), # CLOSING SM BandStructure

         ])
    
    KernelOptimSubMatcher_edft = SM(name= 'energy_components',
                                    startReStr = r"\>\>\> Density kernel optimised for the current NGWF basis\:",
                                    sections = ['x_onetep_section_kernel_optimisation'],
                                    
                                    repeats = True,

                subMatchers = [                     
                    SM(r"\s*Total free energy\s*\=\s*(?P<x_onetep_total_free_energy>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*Total energy\s*\=\s*(?P<x_onetep_total_energy>[-+0-9.eEdD]*)\s*"), # matching final converged total energy
                    SM(r"\s*Estimated bandgap\s*\=\s*(?P<x_onetep_band_gap>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*RMS occupancy error\s*\=\s*(?P<x_onetep_rms_occupancy_error>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*\[H\,K\] commutator\s*\=\s*(?P<x_onetep_commutator>[-+0-9.eEdD]*)\s*"),
                    ])
   
    KernelOptimSubMatcher = SM(name= 'energy_components',
                                    startReStr = r"\>\>\> Density kernel optimised for the current NGWF basis\:",
                                    sections = ['x_onetep_section_kernel_optimisation'],
                                    
                                    repeats = True,

                subMatchers = [                     
                    SM(r"\s*Total free energy\s*\=\s*(?P<x_onetep_total_free_energy>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*Total energy\s*\=\s*(?P<x_onetep_total_energy>[-+0-9.eEdD]*)\s*"), # matching final converged total energy
                    SM(r"\s*Estimated bandgap\s*\=\s*(?P<x_onetep_band_gap>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*RMS occupancy error\s*\=\s*(?P<x_onetep_rms_occupancy_error>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*\[H\,K\] commutator\s*\=\s*(?P<x_onetep_commutator>[-+0-9.eEdD]*)\s*"),
                    ])


    energycomponentsSubMatcher = SM(name= 'energy_components',
                                    startReStr = r"\s*\-\-\-\-\-\-\-\-\-\-* ENERGY COMPONENTS \(Eh\) \-\-\-\-\-\-\-\-\-\-\-*",
                                    # endReStr = r"\sBFGS\:\sfinished iteration\s*\0\s*",
                                    sections = ['x_onetep_section_energy_components'],
                                    forwardMatch = True,

                subMatchers = [                     
                    
                    SM(r"\s*\| Kinetic\s*\:\s*(?P<x_onetep_electronic_kinetic_energy>[-+0-9.eEdD]*)\s\|\s*"), # matching final converged total energy
                    SM(r"\s*\| Pseudopotential \(local\)\s*\:\s*(?P<x_onetep_pseudo_local_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Pseudopotential \(non\-local\)\s*\:\s*(?P<x_onetep_pseudo_non_local_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Hartree\s*\:\s*(?P<x_onetep_energy_correction_hartree_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Exchange\-correlation\s*\:\s*(?P<x_onetep_energy_XC_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Ewald\s*\:\s*(?P<x_onetep_ewald_correction_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Dispersion Correction\s*\:\s*(?P<x_onetep_dispersion_correction_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*Integrated density\s*\:\s*(?P<x_onetep_integrated_density_store>[-+0-9.eEdD]*)\s*"), 

                    ])
    
    edft_SubMatcher = SM (name = 'EDFT submatcher',
            sections = ['x_onetep_section_edft'],
            startReStr = r"\-\-\-\-* Ensemble\-DFT optimisation \-\-\-\-*", 
            repeats = True,
            # endReStr = r"\-*\s*",
            subMatchers = [ 
                SM (name = 'EDFT submatcher',
                sections = ['x_onetep_section_edft_iterations'],
                startReStr = r"\#\s*(?P<x_onetep_edft_iteration>[0-9]+)\s*(?P<x_onetep_edft_rms_gradient>[-\d\.]+)\s*(?P<x_onetep_edft_commutator>[\d\.]+)\s*(?P<x_onetep_edft_free_energy>[\d\.]+)",
                repeats = True,
                # endReStr = r"\-*\s*",
                    subMatchers = [ 
                        SM(r"\#\s*(?P<x_onetep_edft_iteration>[0-9]+)\s*(?P<x_onetep_edft_rms_gradient>[-\d\.]+)\s*(?P<x_onetep_edft_commutator>[\d\.]+)\s*(?P<x_onetep_edft_free_energy>[\d\.]+)"),  
                        SM(r"Step\s*\=\s*(?P<x_onetep_edft_step>[\d\.]+)\s*"), 
                        SM(r"Energy\s*\=\s*(?P<x_onetep_edft_energy>[\d\.]+)\s*"), 
                        SM(r"Est\. 0K Energy 0\.5\*\(E\+A\)\s*\=\s*(?P<x_onetep_edft_0K>[\d\.]+)\s*"), 
                        SM(r"Residual Non\-orthogonality \=\s*(?P<x_onetep_residual_nonorthog>[\d\.]+)\s*"), 
                        SM(r"Residual N\_electrons\s*\=\s*(?P<x_onetep_residual_n_elec>[\d\.]+)\s*"), 
                        SM (name = 'EDFT submatcher',
                            sections = ['x_onetep_section_edft_spin'],
                            startReStr = r"\s*(?P<x_onetep_edft_spin_type>[0-9]+)\s*(?P<x_onetep_edft_n_electrons>[\d\.]+)\s*(?P<x_onetep_edft_fermi_level>[-\d\.]+)\s*(?P<x_onetep_edft_fermi_level_delta>[-\d\.]+)",
                            repeats = True,
                            # endReStr = r"\s*\-\-\-\-\-\-\s*",
                            subMatchers = [ 
                                     SM (name = 'EDFT submatcher',
                                        sections = ['x_onetep_section_edft_spin_iterations'],
                                        startReStr = r"\s*(?P<x_onetep_edft_orbital_iteration_spin>[0-9]+)\s\|\s*(?P<x_onetep_edft_eigenvalue>[-\d\.]+)\s*(?P<x_onetep_edft_occupancy>[-\d\.]+)\s\|",
                                        repeats = True,
                                                
                                            )]),
                        KernelOptimSubMatcher_edft,
                     
                     ]) ])
                
    KernelOptimSubMatcher_frame = SM(name= 'energy_components',
                                    startReStr = r"\>\>\> Density kernel optimised for the current NGWF basis\:",
                                    sections = ['x_onetep_section_kernel_optimisation'],
                                    
                                    repeats = True,

                subMatchers = [                     
                    SM(r"\s*Total free energy\s*\=\s*(?P<x_onetep_total_free_energy>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*Total energy\s*\=\s*(?P<x_onetep_total_energy>[-+0-9.eEdD]*)\s*"), # matching final converged total energy
                    SM(r"\s*Estimated bandgap\s*\=\s*(?P<x_onetep_band_gap>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*RMS occupancy error\s*\=\s*(?P<x_onetep_rms_occupancy_error>[-+0-9.eEdD]*)\s*"),
                    SM(r"\s*\[H\,K\] commutator\s*\=\s*(?P<x_onetep_commutator>[-+0-9.eEdD]*)\s*"),
                    ])

    
    energycomponentsSubMatcher_frame = SM(name= 'energy_components',
                                    startReStr = r"\s*\-\-\-\-\-\-\-\-\-\-* ENERGY COMPONENTS \(Eh\) \-\-\-\-\-\-\-\-\-\-\-*",
                                    # endReStr = r"\sBFGS\:\sfinished iteration\s*\0\s*",
                                    sections = ['x_onetep_section_energy_components'],
                                    forwardMatch = True,

                subMatchers = [                     
                    
                    SM(r"\s*\| Kinetic\s*\:\s*(?P<x_onetep_electronic_kinetic_energy>[-+0-9.eEdD]*)\s\|\s*"), # matching final converged total energy
                    SM(r"\s*\| Pseudopotential \(local\)\s*\:\s*(?P<x_onetep_pseudo_local_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Pseudopotential \(non\-local\)\s*\:\s*(?P<x_onetep_pseudo_non_local_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Hartree\s*\:\s*(?P<x_onetep_energy_correction_hartree_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Exchange\-correlation\s*\:\s*(?P<x_onetep_energy_XC_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Ewald\s*\:\s*(?P<x_onetep_ewald_correction_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*\| Dispersion Correction\s*\:\s*(?P<x_onetep_dispersion_correction_store>[-+0-9.eEdD]*)\s\|\s*"), 
                    SM(r"\s*Integrated density\s*\:\s*(?P<x_onetep_integrated_density_store>[-+0-9.eEdD]*)\s*"), 

                    ])
    
    singlepointSubMatcher = SM(name = 'single_point',
                # startReStr = r"\s*\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\- ENERGY COMPONENTS \(Eh\) \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-" or              
                startReStr = r"\s*\<\<\<\<\< CALCULATION SUMMARY \>\>\>\>\>\s*",
                # startReStr = r"\>\>\> Density kernel optimised for the current NGWF basis\:",
                forwardMatch = True,
                # startReStr = r"\s*\<* CALCULATION SUMMARY \>*\s*", 
                endReStr = r"\sBFGS\:\sfinished iteration\s*0\s",
                sections = ["section_single_configuration_calculation"],
                subMatchers = [ 
                    

                    
                    SM(sections = ['section_scf_iteration'],
                        startReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*\<\-\-\sCG\s*"),
                        # endReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*\<\-\-\sCG\s*",
                        # endReStr = r"\s*[0-9]+\s*[+0-9.eEdD]+\s*[-+0-9.eEdD]*\s*\<\-\-\sCG\s*",
                        # subMatchers = [   
                        #     SM(r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*\<\-\-\sCG\s*",
                        #     repeats = True)]),
                          
                    SM(sections = ['section_scf_iteration'],
                        startReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*(?P<energy_change_scf_iteration>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*",repeats = True,
                        # endReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*\<\-\-\sCG\s*",
                        endReStr = r"\s*[0-9]+\s*[+0-9.eEdD]+\s*[-+0-9.eEdD]*\s*\<\-\-\sCG\s*",
                        # subMatchers = [   
                        #     SM(r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*\<\-\-\sCG\s*",
                        #     repeats = True)]),
                          ),  
                    SM(r"\<QC\>\s*\[NGWF iterations]\:\s*(?P<x_onetep_n_ngwf_iterations>[0-9]*)\s*"),
                    SM(r"\<QC\>\s*\[total\_energy\]\:\s*(?P<energy_total>[-+0-9.eEdD]*)\s*"), # matching final converged total energy
                    SM(r"\<QC\>\s*\[rms\_gradient\]\:\s*(?P<x_onetep_final_rms_gradient>[-+0-9.eEdD]*)\s*"),
    
                 
                    SM(startReStr = r"\*\*\*\*\** Forces \*\*\*\*\**\s*",
                         subMatchers = [
                                    
                                    SM(r"\*\s*[A-Za-z]+\s*[0-9]+\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\s\*",
                                        repeats = True)
                         ]), 
                    
                    # SM(#startReStr = r"\sBFGS\:\sFinal\sbulk\smodulus\sunchanged\sfrom\sinitial\svalue\s*",
                    #     startReStr = r"\s\*\*\*\*\**\sSymmetrised Forces\s\*\*\*\*\**\s*",
                    #     subMatchers = [
                    #        SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                    #           repeats = True)
                    #                   ]),
                     
                    # SM(name = 'stresstensor',
                    #     #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                    #     startReStr = r"\s\*\*\*\*\** Stress Tensor \*\*\*\*\**\s*",
                    #     #GGGrepeats = True,
                    #     sections = ['x_onetep_section_stress_tensor'],
                    #     subMatchers = [
                    #        SM(r"\s*\*\s*[a-z]\s*(?P<x_onetep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                    #           repeats = True),  
                    #                   ]), # CLOSING section_stress_tensor
                   
                     
                    # SM(name = 'stresstensor',
                    #     #startReStr = r"\s\*\*\*\*\*\*\*\*\*\ Symmetrised Stress Tensor \*\*\*\*\*\*\*\*\*\*\*\s*",
                    #     startReStr = r"\s\*\*\*\*\** Symmetrised Stress Tensor \*\*\*\*\**\s*",
                    #     sections = ['x_onetep_section_stress_tensor'],
                    #     subMatchers = [
                    #        SM(r"\s*\*\s*[a-z]\s*(?P<x_onetep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                    #           repeats = True),  
                    #                   ]), # CLOSING section_stress_tensor
                    # SM(startReStr = r"\s*x\s*MD\sData\:\s*x",
                    #      subMatchers = [
                    #         SM(r"\s*x\s*time\s*\:\s*(?P<x_onetep_frame_time_0>[+0-9.eEdD]+)\s*ps\s*x\s*"),
                    # ]),

                 ])     
    
    MDSubMatcher = SM(name = 'MD',
                startReStr = r"\sStarting MD iteration\s*[0-9.]+\s\.\.\.\s*",             
                #endReStr = r"BFGS\: finished iteration     0 with enthalpy\= \-2\.14201700E\+002 eV",
                sections = ["x_onetep_section_SCF_iteration_frame"],
                # sections = ["section_single_configuration_calculation"],
                endReStr = r"\s\.\.\.\sfinished MD iteration\s*[0-9.]+\s*",
                repeats =True,
                subMatchers = [                     
                    # SM(r"\s*[0-9]+\s*(?P<onetep_SCF_frame_energy>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                    #     endReStr = "\n",
                    #     repeats = True),
                    SM(r"\s*[0-9]+\s*(?P<x_onetep_SCF_frame_energy>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*(?P<x_onetep_SCF_frame_energy_gain>[-+0-9.eEdD]*)\s*(?P<x_onetep_frame_time_scf_iteration_wall_end>[0-9.]*)\s*\<\-\-\sSCF\s*",
                        endReStr = "\n",
                        repeats = True),                
                    SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<x_onetep_frame_energy_free>[-+0-9.eEdD]*)"), # matching final converged total free energy
                    SM(r"NB est\. 0K energy\s*\(E\-0\.5TS\)\s*= *(?P<x_onetep_frame_energy_total_T0>[-+0-9.eEdD]*)"), # 0K corrected final SCF energy
                    SM(startReStr = r"\s*x\s*MD\sData\:\s*x",
                         subMatchers = [
                            SM(r"\s*x\s*time\s*\:\s*(?P<x_onetep_frame_time>[+0-9.eEdD]+)\s*ps\s*x\s*"),
                         ]),
                   

                    ])
    ########################################
    # Sub matcher for geometry optimisation 
    ########################################

    geomOptimSubMatcher_improving_cycle =  SM (name = 'geometry_optimisation_improving',
            startReStr = r"\s[A-Za-z]+\:\simproving\siteration\s*1\s+",#with\sline\sminimization\s\(lambda\=\s*1\.557176\)\s*",
       
            endReStr = r"\s[A-Za-z]+\:\sfinished iteration\s*1",
       
            subMatchers = [               
              
                        SM(name = 'cellInformation',
                            startReStr = r"\s*Unit Cell\s*",
                            forwardMatch = True,
                            sections = ["x_onetep_section_cell_optim"],
                            subMatchers = [
                                SM(r"\s*(?P<x_onetep_cell_vector_optim>[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEd]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                                    endReStr = "\n",
                                    repeats = True),

                             ]), # CLOSING onetep_section_cell


           # atomic positions and cell dimesions
                        SM(startReStr = r"\s*Lattice parameters",
                            forwardMatch = True,
                            sections = ["x_onetep_section_atom_positions_optim"],
                            subMatchers = [

                                SM(r"\s*a \=\s*(?P<x_onetep_cell_length_a_optim>[\d\.]+)\s*alpha \=\s*(?P<x_onetep_cell_angle_alpha_optim>[\d\.]+)"),
                                SM(r"\s*b \=\s*(?P<x_onetep_cell_length_b_optim>[\d\.]+)\s*beta  \=\s*(?P<x_onetep_cell_angle_beta_optim>[\d\.]+)"),
                                SM(r"\s*c \=\s*(?P<x_onetep_cell_length_c_optim>[\d\.]+)\s*gamma \=\s*(?P<x_onetep_cell_angle_gamma_optim>[\d\.]+)"),

                            ]), # CLOSING onetep_section_atom_positions

                        SM(r"\s*x\s*(?P<x_onetep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<x_onetep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                            endReStr = "\n",
                            repeats = True),

                        
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*(?P<time_scf_iteration_wall_end>[0-9.]*)\s*\<\-\-\sSCF\s*",
                             endReStr = "\n",
                             repeats = True),

                                   # ]), # CLOSING section_scf_iteration                                      
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*(?P<time_scf_iteration_wall_end>[0-9.]*)\s*\<\-\-\sSCF\s*",
                            endReStr = "\n",
                            repeats = True),           
                            
                        SM(r"Final energy = *(?P<x_onetep_improved_energy_total>[-+0-9.eEdD]*)"), # matching final converged total energy
                        SM(r"Final energy\,\s*E\s*= *(?P<x_onetep_improved_energy_total>[-+0-9.eEdD]*)"), # matching final converged total energy
                 #SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<onetep_energy_free>[-+0-9.eEdD]*)"),
                

                        SM(name = 'stresstensor',
                         #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                        startReStr = r"\s\*\*\*\*\** Stress Tensor \*\*\*\*\**\s*",
                        forwardMatch = True,
                        #repeats = True,
                        sections = ['x_onetep_section_stress_tensor'],
                        subMatchers = [
                            SM(r"\s*\*\s*[a-z]\s*(?P<x_onetep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                               repeats = True),  
                                       ]), # CLOSING section_stress_tensor

                        
                

                 ])
    
    geomOptimSubMatcher =  SM (name = 'geometry_optimisation',
            startReStr = r"\sStarting BFGS iteration\s*(?P<x_onetep_geom_iteration_index>[0-9.]+)\s\.\.\.\s*",
            sections = ['section_single_configuration_calculation','section_system'],
            endReStr = r"\sBFGS\s\:\sGeometry optimization [a-z]+ to converge after\s*",
            #endReStr = r"\s\[A-Za-z]+\:\sGeometry\soptimization\scompleted\ssuccessfully.\s*",
            repeats = True,
            subMatchers = [               


                        SM(r"\s*x\s*(?P<x_onetep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<x_onetep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                            endReStr = "\n",
                            repeats = True),

                        KernelOptimSubMatcher_frame,
                        
                        energycomponentsSubMatcher_frame,
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_rms_gradient>[+0-9.eEdD]+)\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)\s*(?P<energy_change_scf_iteration>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*",repeats = True,
                            endReStr = r"\s*[0-9]+\s*[+0-9.eEdD]+\s*[-+0-9.eEdD]*\s*\<\-\-\sCG\s*"),  
                        
                        SM(r"\<QC\>\s*\[NGWF iterations]\:\s*(?P<x_onetep_n_ngwf_iterations>[0-9]*)\s*"),
                        SM(r"\<QC\>\s*\[total\_energy\]\:\s*(?P<energy_total>[-+0-9.eEdD]*)\s*"), # matching final converged total energy
                        SM(r"\<QC\>\s*\[rms\_gradient\]\:\s*(?P<x_onetep_final_rms_gradient>[-+0-9.eEdD]*)\s*"),
    
                    
                        SM(startReStr = r"\*\*\*\*\** Forces \*\*\*\*\**\s*",
                            subMatchers = [
                                    SM(r"\*\s*[A-Za-z]+\s*[0-9]+\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\s\*",
                                        repeats = True)
                         ]), 

                        # geomOptimSubMatcher_improving_cycle,
                        SM(r"\sBFGS\s\:\sGeometry optimization (?P<x_onetep_geom_converged>[a-z]+) to converge after\s*"),
                        SM(r"\s[A-Za-z]+\:\sGeometry\soptimization\scompleted\s(?P<x_onetep_geom_converged>[a-z]+)\.\s*"),
                        ])
    
    
    
    geomOptim_finalSubMatcher = SM (name = 'geometry_optimisation_final_configuration',
            startReStr = r"\s[A-Za-z]+\:\sFinal Configuration\:",
            sections = ['section_single_configuration_calculation','section_system'],

            repeats = True,
            subMatchers = [                  
                                
                        

                            SM(r"\s*x\s*(?P<x_onetep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<x_onetep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                                    endReStr = "\n",
                                    repeats = True),
                            

                        
                            # SM(r"\s[A-Za-z]+\:\sFinal\sEnthalpy\s*\=\s(?P<x_onetep_enthalpy__eV>[-+0-9.eEdD]+)"),        
                            # SM(r"\s[A-Za-z]+\:\sFinal\s\<frequency\>\s*\=\s*(?P<x_onetep_frequency>[-+0-9.eEdD]+)"),                    
                            # SM(#startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* Forces \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
                            #         startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                                    
                            #         subMatchers = [
                            #            SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                            #               repeats = True)
                            #                       ]),
                            # SM(#startReStr = r"\sBFGS\:\sFinal\sbulk\smodulus\sunchanged\sfrom\sinitial\svalue\s*",
                            #         startReStr = r"\s\*\*\*\*\**\sSymmetrised\sForces\s\*\*\*\*\**\s*",
                            #         repeats = True,
                            #         subMatchers = [
                            #            SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                            #               repeats = True)
                            #                       ]), 
                            
                            # SM(name = 'Forces',
                            #     startReStr = r"\s\*\*\*\*\** TDDFT Forces \*\*\*\*\**\s*",
                            #     subMatchers = [
                            #         SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                            #             repeats = True)
                            #               ]),

                            # SM(name = 'Forces',
                            #     startReStr = r"\s\*\*\*\*\** TDDFT Symmetrised Forces \*\*\*\*\**\s*",
                            #     subMatchers = [
                            #         SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<x_onetep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                            #             repeats = True)
                            #               ]),
                            
                            
                   

                              ])
                
    
    TSSubMatcher =  SM (name = 'TS_search',
            startReStr = r"\+* Starting Transition State Search \+*",
            sections = ['section_single_configuration_calculation','section_system'],
            endReStr = r"\s*Atomic\sPopulations\s\(Mulliken\)\s*",
            
            repeats = True,
            subMatchers = [
                        SM (name = 'TS_search',
                        startReStr = r"SCF loop\s*Energy\s*Fermi\s*Energy gain\s*Timer\s*\<\-\-\sSCF\s*",
                        sections = ['x_onetep_section_ts_scf'],
                        endReStr = r"\sPath coordinate\:",
            #endReStr = r"\s\[A-Za-z]+\:\sGeometry\soptimization\scompleted\ssuccessfully.\s*",
                        repeats = True,
                        subMatchers = [
                            SM(sections = ['x_onetep_section_ts_scf_iteration'],
                                startReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_ts_iteration_energy__eV>[-+0-9.eEdD]*)\s*(?P<x_onetep_scf_ts_iteration_energy_change__eV>[-+0-9.eEdD]*)\s*(?P<x_onetep_scf_ts_time>[0-9.]*)\s*\<\-\-\sSCF\s*",
                                endReStr = "\n",
                                repeats = True),
                        
                            SM(sections = ['x_onetep_section_ts_scf_iteration'],
                                startReStr = r"\s*[0-9]+\s*(?P<x_onetep_scf_ts_iteration_energy__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*(?P<x_onetep_scf_ts_iteration_energy_change__eV>[-+0-9.eEdD]*)\s*(?P<x_onetep_scf_ts_time>[0-9.]*)\s*\<\-\-\sSCF\s*",
                                endReStr = "\n",
                                repeats = True),     
                            
                        SM(r"Final energy = *(?P<x_onetep_scf_ts_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                        SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<x_onetep_scf_ts_total_energy_free__eV>[-+0-9.eEdD]*)"),
                        SM(r"NB est\. 0K energy\s*\(E\-0\.5TS\)\s*= *(?P<x_onetep_scf_ts_T0__eV>[-+0-9.eEdD]*)"),
                
                 ]),   
             ])
    
    Mulliken_SubMatcher = SM (name = 'Mulliken population analysis',
            startReStr = r"\s*Mulliken Atomic Populations\s*",
            sections = ['x_onetep_section_population_analysis'],
            endReStr = r"\s*Bond\s*Population\s*Spin\s*Length\s\(bohr\)\s*",
            #endReStr = r"\s*Bond\s*Population\s*Length\s\(A\)\s*",
            repeats = True,
            subMatchers = [ 
                SM(r"\s*[a-zA-Z]+\s*[0-9.]+\s*(?P<x_onetep_total_orbital_store>[0-9.]+)\s*(?P<x_onetep_mulliken_charge_store>[0-9.]+)\s*", repeats = True),
                SM(r"\s*[a-zA-Z]+\s*[0-9.]+\s*(?P<x_onetep_total_orbital_store>[0-9.]+)\s*(?P<x_onetep_mulliken_charge_store>[0-9.]+)\s*(?P<x_onetep_spin_store>[0-9.]+)\s*",       
                    
                    repeats = True),
                 ]) 
    Orbital_SubMatcher = SM (name = 'Orbital Information',
            startReStr = r"\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\= Orbital energy and occupancy information \=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\s*",
            sections = ['x_onetep_section_orbital_information'],
            forwardMatch = True,
            # endReStr = "\n",
            endReStr =r"\s*\.\.\.\.\.\.\.\s*\-\-\- gap \-\-\s*\.\.\.\.\.\.\.\.\.\s*",
            
            subMatchers = [ 
                SM(r"\s*Total number of orbitals\:\s*(?P<x_onetep_total_number_orbitals>[0-9]+)\s*"),
                SM(r"\s*Number of occupied orbitals\:\s*(?P<x_onetep_total_number_occ_orbitals>[0-9]+)\s*"),
                SM(r"\s*Occupancy sum\:\s*(?P<x_onetep_occupancy_sum>[\d\.]+)\s*"),
                SM(r"\s*HOMO\-LUMO gap\:\s*(?P<x_onetep_homo_lumo_gap>[\d\.]+)\s*"),
                SM(r"\s*Mid\-gap level\:\s*(?P<x_onetep_mid_gap>[\d\.]+)\s*"), 
                SM(r"\s*(?P<x_onetep_orbital_number>[0-9]+)\s*(?P<x_onetep_orbital_energy>[-\d\.]+)\s*(?P<x_onetep_orbital_occupancy>[\d\.]+)\s*",repeats=True),  
                 ]) 
    Orbital_SubMatcher_2 = SM (name = 'Orbital Information',
            startReStr = r"\s*\.\.\.\.\.\.\.\s*\-\-\- gap \-\-\s*\.\.\.\.\.\.\.\.\.\s*",
            sections = ['x_onetep_section_orbital_information'],
            forwardMatch = True,
            endReStr = r"\s*\.\.\.\.\.\.\.\s*\.*\s*\.\.\.\.\.\.\.\.\.\s*",
            # endReStr =r"\s*\.\.\.\.\.\.\.\s*\-\-\- gap \-\-\s*\.\.\.\.\.\.\.\.\.\s*",
            
            subMatchers = [ 
                SM(r"\s*(?P<x_onetep_orbital_number>[0-9]+)\s*(?P<x_onetep_orbital_energy>[-\d\.]+)\s*(?P<x_onetep_orbital_occupancy>[\d\.]+)\s*",repeats=True),  
                 ]) 
    Orbital_SubMatcher_3 = SM (name = 'Orbital Information',
            startReStr = r"\s*\.\.\.\.\.\.\.\s*\.*\s*\.\.\.\.\.\.\.\.\.\s*",
            sections = ['x_onetep_section_orbital_information'],
            
            endReStr = "\n",
            # endReStr =r"\s*\.\.\.\.\.\.\.\s*\-\-\- gap \-\-\s*\.\.\.\.\.\.\.\.\.\s*",
            
            subMatchers = [ 
                SM(r"\s*(?P<x_onetep_orbital_number>[0-9]+)\s*(?P<x_onetep_orbital_energy>[-\d\.]+)\s*(?P<x_onetep_orbital_occupancy>[\d\.]+)\s*",repeats=True),  
                 ]) 
    

    ########################################
    # return main Parser ###################
    ########################################
    return SM (name = 'Root',
        startReStr = "",
        forwardMatch = True,
        # subFlags = SM.SubFlags.Unordered,
        weak = True,
        subMatchers = [
        SM (name = 'NewRun',
            startReStr = r"\s\|\s*Linear-Scaling Ab Initio Total Energy Program\s*\|\s*",
            subFlags = SM.SubFlags.Unordered,
            required = True,
            forwardMatch = True,
            sections = ['section_run'],
            subMatchers = [
                
    
            
                SM(name = 'ProgramHeader',
                  startReStr = r"\s\|\s*Linear-Scaling Ab Initio Total Energy Program\s*\|\s*",
                  subMatchers = [

                     SM(r"\s\|\s*Version\s(?P<program_version>[0-9a-zA-Z_.]*)"),
                     SM(r"\s\|\s*in all publications arising from your use of (?P<program_name>[a-zA-Z]+)*"),
                     SM(r"\s*Default threads\:\s(?P<x_onetep_number_of_processors>[0-9.]*)"),
                     # SM(r"\sCompiler\: *(?P<x_onetep_compiler>[a-zA-Z\s0-9.]*)"),
                     # SM(r"\sMATHLIBS\: *(?P<x_onetep_maths_library>[a-zA-Z0-9.() ]*)\s*"),
                     # SM(r"\sFFT Lib \: *(?P<x_onetep_fft_library>[a-zA-Z0-9.() ]*)\s*"),
                     # SM(r"\sFundamental constants values\: *(?P<x_onetep_constants_reference>[a-zA-Z0-9.() ]*)\s*"),

                                  ]), # CLOSING SM ProgramHeader
                                
               
                # scfEigenvaluesSubMatcher, # export section_eigenvalues_group to the correct nesting


                 # section_method

                calculationMethodSubMatcher,
                basisSetCellAssociatedSubMatcher,
      
                SM(name = 'Atom_topology',
                  startReStr = r"\%block\s*species\s*",              
                  
                  # endReStr = r"\%endblock species\s*",  
                  #forwardMatch = True,
                  sections = ['section_topology'],
                  subMatchers = [
                 
                      SM(r"[0-9a-zA-Z]+\s*(?P<x_onetep_store_atom_name>[a-zA-Z]+)\s*(?P<x_onetep_store_atom_mass>[0-9.]+)\s*(?P<x_onetep_n_ngwf_atom_store>[0-9.]+)\s*(?P<x_onetep_ngwf_radius_store>[\d\.]+)\s*",
                        #endReStr = "\n",   
                        repeats = True),
                     ]), # CLOSI
                ElectronicParameterSubMatcher,
                
                
                # ElectronicMinimisParameterSubMatcher,

                # DensityMixingParameterSubMatcher,

                # PopulationAnalysisParameterSubMatcher,

                # MDParameterSubMatcher,
                
                # GeomOptimParameterSubMatcher,
                
                # phononCalculationSubMatcher,
                
                # TSParameterSubMatcher,
                
                # OpticsParameterSubMatcher,

                # ElecSpecParameterSubMatcher,
                
                # TDDFTParameterSubMatcher,

                # systemDescriptionSubMatcher, # section_system subMatcher
        
                
               

          #       SM(r"\s*Point group of crystal\s\=\s*[0-9.]+\:\s(?P<x_onetep_crystal_point_group>[a-zA-Z0-9.]+)"),
          #       SM(r"\s*Space group of crystal\s\=\s*[0-9.]+\:\s(?P<x_onetep_space_group>[a-zA-Z0-9.]+)"),

          # ############ onetep-specific van der Waals method parameters #############################     
                
          #       SM(name = "van der Waals onetep TS",
          #          startReStr = r"\s*Dispersion\-correction scheme\s\:\s+",
          #          forwardMatch = True,
          #          sections = ["x_onetep_section_van_der_Waals_parameters"],
          #          subMatchers = [
          #               ######## Method TS #######
          #               SM(r"\s*Parameter sR\s*\: *(?P<x_onetep_Parameter_sR> [0-9.]+)"),
          #               SM(r"\s*Parameter d\s*\: *(?P<x_onetep_Parameter_d> [0-9.]+)"),
          #               ######## Method OBS #######
          #               SM(r"\s*Parameter lambda\s*\: *(?P<x_onetep_Parameter_LAMBDA> [0-9.]+)"),
          #               SM(r"\s*Parameter n\s*\: *(?P<x_onetep_Parameter_n> [0-9.]+)"),
               
          #               ######## Method G06 #######
          #               SM(r"\s*Parameter s6\s*\: *(?P<x_onetep_Parameter_s6> [0-9.]+)"),
          #               SM(r"\s*Parameter d\s*\: *(?P<x_onetep_Parameter_d> [0-9.]+)")
          #                     ]), # CLOSING van der waals onetep parameters
              
          #       #geomOptimSubMatcher_init,
             
                
                                     
          #       SM(r"Calculating total energy with cut\-off of  (?P<x_onetep_basis_set_planewave_cutoff_iteration_0>[0-9.]+)",
          #                      repeats = True,),
                # KernelOptimSubMatcher,
                # energycomponentsSubMatcher,
                edft_SubMatcher,
                KernelOptimSubMatcher,
                energycomponentsSubMatcher,
                singlepointSubMatcher,
                geomOptimSubMatcher,
                # TSSubMatcher,
                # MDSubMatcher,

                # SM(name = "Vibrational_frequencies",
                #     sections = ["x_onetep_section_vibrational_frequencies"],
                #     startReStr = r"\s\+\s*q\-pt\=\s*[0-9.]+",
                #     #startReStr = r"\s\+\s*Vibrational Frequencies\s*\+\s*",
                #     repeats = True,
                #     subMatchers = [
                #         SM(r"\s\+\s*(?P<x_onetep_n_iterations_phonons>[0-9.]+)\s*(?P<x_onetep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<x_onetep_ir_store>[a-z])\s*(?P<x_onetep_ir_intensity_store>[-\d\.]+)\s*[A-Z]\s*(?P<x_onetep_raman_activity_store>[-\d\.]+)\s*[A-Z]\s*\+\s*",       
                #             repeats = True,
                #             endReStr = r"\s\+\s\.*\s\+\s*",
                #              ),
                #         SM(r"\s\+\s*(?P<x_onetep_n_iterations_phonons>[0-9.]+)\s*(?P<x_onetep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<x_onetep_ir_store>[a-z])\s*(?P<x_onetep_ir_intensity_store>[-\d\.]+)\s*[A-Z]\s*[A-Z]\s*\+\s*",       
                #             repeats = True,
                #             endReStr = r"\s\+\s\.*\s\+\s*",
                #              ),
                #         SM(r"\s\+\s*(?P<x_onetep_n_iterations_phonons>[0-9.]+)\s*(?P<x_onetep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<x_onetep_ir_store>[a-z])\s*\+\s*",          
                #             repeats = True,
                #             endReStr = r"\s\+\s\.*\s\+\s*",
                #             ),
                #           ]),
                
                # bandStructureSubMatcher,  # band structure subMatcher
                                       
                # SM(name = 'ramantensor',
                #         #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                #         startReStr = r"\s\+\sMode number\:\s*[0-9.]+\sRaman\stensor\s*Depolarisation Ratio\s*\+\s*",
                #         endReStr = r"\s\+\s*\+\s*",
                #         repeats = True,
                #         sections = ['x_onetep_section_raman_tensor'],
                #         subMatchers = [
                #            SM(r"\s\+\s*(?P<x_onetep_store_raman_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                #               repeats = True),  
                #                       ]), # CLOSIN

                # SM(startReStr = r"\s\*\*\*\*\**\sSymmetrised Forces\s\*\*\*\*\**\s*",
                #         subMatchers = [
                #            SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<onetep_store_atom_forces_band>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                #               repeats = True)
                #                       ]),
                # SM(startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                #     subMatchers = [
                #            SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<onetep_store_atom_forces_band>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                #               repeats = True)
                #                       ]),

                     
                
                    
                # geomOptimSubMatcher,
                # # geomOptimSubMatcherDMD,
                # # geomOptimSubMatcherDI,
                geomOptim_finalSubMatcher,            
                
                Mulliken_SubMatcher,
                
                Orbital_SubMatcher,
                
                Orbital_SubMatcher_2 ,
                
                Orbital_SubMatcher_3,
                SM(name = 'calc_time',
                    startReStr = r"\-\-\-*\sTIMING INFORMATION\s\-\-\-*",
                    sections = ['x_onetep_section_time'],
                    subMatchers = [
                        SM(r"AVERAGE TIME\:\s*(?P<x_onetep_avarage_time>[0-9.]*)"), # matching final converged total energy
                        SM(r"TOTAL TIME\:\s*(?P<x_onetep_total_time>[0-9.]*)"),
                        SM(r"Job completed\:\s*(?P<x_onetep_final_date>[0-9.-]*)(?P<x_onetep_final_time>[0-9.:]*)"),
                                      ]), # CLOSING section_onetep_time
                                        
                  
                
                         ])       
                
 
        ]) # CLOSING SM NewRun 



def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = { 'x_onetep_ir_intensity_store' : CachingLevel.Cache,
                                'x_onetep_ir_store' : CachingLevel.Cache,
                                'x_onetep_vibrationl_frequencies_store' : CachingLevel.Cache,
                                'x_onetep_n_iterations_phonons' : CachingLevel.Cache,
                                'x_onetep_mulliken_charge_store' : CachingLevel.Cache,
                                'x_onetep_orbital_contributions' : CachingLevel.Cache,
                                # 'x_onetep_section_cell_optim': CachingLevel.Cache,
                                'x_onetep_section_atom_positions_optim' : CachingLevel.Cache,
                                'x_onetep_section_eigenvalues':CachingLevel.Cache,
                                'x_onetep_section_k_points':CachingLevel.Cache,
                                'x_onetep_section_k_band':CachingLevel.Cache,
                                # 'band_energies' : CachingLevel.Cache,
                                # 'band_k_points' : CachingLevel.Cache,
                                'x_onetep_basis_set_planewave_cutoff' : CachingLevel.Cache,
                                # 'eigenvalues_values': CachingLevel.Cache,
                                # 'eigenvalues_kpoints':CachingLevel.Cache,
                                'x_onetep_total_energy_corrected_for_finite_basis_store': CachingLevel.Cache,
                                'x_onetep_frame_time':CachingLevel.Cache,
                                'x_onetep_section_SCF_iteration_frame':CachingLevel.Cache,
                                'x_onetep_raman_activity_store': CachingLevel.Cache,
                                'x_onetep_SCF_frame_energy_gain':CachingLevel.Cache,
                                'x_onetep_frame_energy_free':CachingLevel.Cache,
                                'x_onetep_frame_energy_total_T0':CachingLevel.Cache,
                                'x_onetep_frame_time_scf_iteration_wall_end':CachingLevel.Cache,
                                'x_onetep_total_orbital_store':CachingLevel.Cache,
                                'x_onetep_mulliken_charge_store':CachingLevel.Cache,
                                'x_onetep_spin_store':CachingLevel.Cache,
                                'x_onetep_ngwf_radius_store':CachingLevel.Cache,
                                'x_onetep_n_ngwf_atom_store':CachingLevel.Cache,
                                'x_onetep_section_energy_components':CachingLevel.Cache,
                                'x_onetep_SCF_frame_energy':CachingLevel.Cache,
                                'x_onetep_electronic_kinetic_energy':CachingLevel.Cache,
                                'x_onetep_pseudo_local_store' :CachingLevel.Cache,
                                'x_onetep_pseudo_non_local_store':CachingLevel.Cache,
                                
                                'x_onetep_energy_correction_hartree_store':CachingLevel.Cache,
                                'x_onetep_energy_XC_store':CachingLevel.Cache,
                                'x_onetep_ewald_correction_store':CachingLevel.Cache,
                                'x_onetep_dispersion_correction_store' :CachingLevel.Cache,
                                'x_onetep_integrated_density_store':CachingLevel.Cache,
                                }

    # Set caching for temparary storage variables
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('x_onetep_store_')
            or name.startswith('x_onetep_cell_')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main():
    """Main function.

    Set up everything for the parsing of the onetep main file and run the parsing.
    """
    # get main file description
    onetepMainFileSimpleMatcher = build_onetepMainFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/onetep.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/onetep.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'onetep-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
    # start parsing
    mainFunction(mainFileDescription = onetepMainFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = OnetepParserContext(),
                 defaultSectionCachingLevel = True)

if __name__ == "__main__":
    main()








