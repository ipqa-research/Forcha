module characterization
   !! This module provides functions to characterize a fluid using input data.
   !! The characterization includes:
   !! - Molecular weight computation
   !! - Density estimation using various methods
   !! - Lumping of components for compositional modeling
   !! - Calculation of critical parameters (Tc, Pc, ω)
   !! The module relies on external modules for data handling 
   !!and property calculations.

   use constants                     
   !! Module containing physical constants
   use dtypes, only: FluidData, FluidDataOut 
   !! Data structures for fluid properties
   use data_from_input, only: data_from_file, FluidData 
   !! Function to read fluid data from a file
   use molecular_weigth              
   !! Module for molecular weight calculations
   use density                        
   !! Module for density calculations
   use lumping                        
   !! Module for lumping methodologies
   use critical_parameters            
   !! Module for critical parameter estimation

contains

   type(FluidDataOut) function characterize(file, mw_source, method, & 
      rho_method, fix_C, eos) result(characterization)
      !! This function performs the full characterization of a reservoir 
      !! fluid based on input data. The characterization process includes 
      !! molecular weight assignment, density computation, lumping of 
      !! components, and estimation of critical parameters.
      !! The function reads fluid composition from a file and applies different
      !! calculation methodologies depending on the specified input parameters.

     !input variables
      character(len=*), intent(in) :: file
      !! Name of the input file containing the fluid composition data.
      character(len=*), intent(in) :: mw_source
      !! Determines whether molecular weights from predefined SCNs in the input 
      !! file are used or whether correlations are used to estimate these values.
      !! Indicates data source ("experimental" or "calculated")
      character(len=*), intent(in), optional :: method
      !! Specifies the method for molecular weight calculation.
      !! Available options:
      !! - "plus_mw": reproduces the molecular weight of the residual fraction.
      !! - "global_mw": reproduces the global molecular weight of fluid.
      character(len=*), intent(in) :: eos
      !! Specifies the equation of state (EOS) used for 
      !! calculating critical parameters.
      logical, intent(in), optional :: fix_C
      !! If present and set to .TRUE., it fixes the value of the constant  'C' 
      !! instead of estimating it dynamically.
      integer, intent(in) :: rho_method
      !! Selector for the density computation method.
      !! Different methods can be used for different fluid 
      !! characterization approaches.
      !
      !! method "1":use experimental scn's densities reported to calculate 
      !! experimental volume for 6plus. The densities values are calculated 
      !! since CN plus
      !
      !! method "2": This densities doesn't experimentals, they can be 
      !! reference values. the densities values are calculated since 6 plus
      !! by correlations.
      !
      !! method "3": use experimental value of density of C6+ to calculate 
      !! the experimental volume
  
      !  Internal variables
      type(FluidData) :: fluid
      !! Structure to store the fluid composition and properties read 
      !! from the input file.
      
      !  Read input fluid data
      fluid = data_from_file(file) 
      !! Read the fluid composition from the input file
      characterization%input_data = fluid 
      !! Store input data in the output structure

      ! Dynamic memory allocation
      allocate(characterization%scn_z(fluid%scn_nc))            
      ! Mole fractions of single carbon numbers (SCNs)
      allocate(characterization%log_scn_z(fluid%scn_nc))        
      ! Logarithm of SCN mole fractions
      allocate(characterization%scn_mw(fluid%scn_nc))           
      ! Molecular weights of SCNs
      allocate(characterization%carbon_number_plus(0))          
      ! Carbon number for the heavy fraction
      allocate(characterization%plus_z_i(0))                    
      ! Mole fraction of the heavy fraction
      allocate(characterization%product_z_mw_plus_i(0))         
      ! Product of mole fraction and molecular weight for the heavy fraction
      allocate(characterization%scn_zm(fluid%scn_nc))           
      ! Mole fraction adjusted for molecular weight calculations
      allocate(characterization%scn_i(0))                       
      ! Indexing for SCN components
      allocate(characterization%plus6_density(0))               
      ! Density of the heavy fraction
      allocate(characterization%mol_fraction(0))                
      ! Mole fractions after characterization
      allocate(characterization%lumped_z(0))                    
      ! Lumped mole fractions
      allocate(characterization%lumped_mw(0))                   
      ! Lumped molecular weights
      allocate(characterization%lumped_densities(0))            
      ! Lumped component densities
      allocate(characterization%critical_temperature(0))        
      ! Estimated critical temperature (Tc)
      allocate(characterization%critical_pressure(0))           
      ! Estimated critical pressure (Pc)
      allocate(characterization%acentric_factor(0))             
      ! Acentric factor (ω)
      allocate(characterization%m_funtion(0))                   
      ! Function for molecular weight adjustment


      ! Computation steps
      ! Step 1: Compute molecular weight for the heavy fraction
      call get_c_or_m_plus(fluid=fluid, mw_source=mw_source, method=method, &
         fix_C=fix_C, characterization=characterization)

      ! Step 2: Compute density using the selected method
      call density_funtion(fluid=fluid, rho_method=rho_method, &
         characterization=characterization)

      ! Step 3: Apply lumping methodology to reduce component count
      call lump(fluid=fluid, characterization=characterization)

      !! Step 4: Compute critical parameters (Tc, Pc, ω) based on EOS
      call get_critical_constants( characterization=characterization, eos=eos)

   end function characterize

end module characterization
