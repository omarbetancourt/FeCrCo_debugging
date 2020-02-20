[Mesh]
  type = GeneratedMesh
  dim = 2
  elem_type = QUAD4
  nx = 80
  ny = 80 # <- may need to reduced if running on local machine
  xmax = 160 # Koyama 2004, Table 1: "Calculation area, LxL = 160nm x 160nm"
  ymax = 160
[]

[Variables] #c1 = Fe, c2 = Cr, c3 = Co

  [c2]  # Mole fraction of Cr (unitless)
  []

  [w2]  # Chemical potential of Cr (eV/mol)
  []

  [c3]  # Mole fraction of Co (unitless)
  []

  [w3]  # Chemical potential of Co (eV/mol)
  []
[]

[AuxVariables]
  [c1]
  []
[]

[AuxKernels]
  [calculate_c1]
    type = ParsedAux
    variable = c1
    function = '1.00 - (c2 + c3)' # Function for computing c1 (composition of Fe)
    args = 'c2 c3'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[ICs] # Initial Conditions for the concentration variables
# Koyama 2004, Fig 3 & 4 = Fe-40at.%Cr-40at.%Co
  [c_CrIC]
    type = RandomIC
    variable = c2
    min = 0.398
    max = 0.402
  []

  [c_CoIC]
    type = RandomIC
    variable = c3
    min = 0.398
    max = 0.402
  []
[]

[BCs] # Boundary Conditions
  [Periodic]
    [c_bcs]
      auto_direction = 'x y'
    []
  []
[]

[Kernels]
# Implementing Off-diagonal Onsager Matrix with Koyama 2004, Equation 8:
# (Example shown in: https://mooseframework.org/source/kernels/SplitCHWRes.html)
  [c2_res]
    type = ADSplitCHParsed
    variable = c2
    f_name = F_total # G_sys
    kappa_name = kappa_c
    w = w2
    args = 'c3 c1'
  []

  [w22_res]
    type = ADSplitCHWRes
    variable = w2
    mob_name = L_22
  []

  [w23_res]
    type = ADSplitCHWRes
    variable = w2
    w = w3
    mob_name = L_23
  []

  [c3_res]
    type = ADSplitCHParsed
    variable = c3
    f_name = F_total
    kappa_name = kappa_c
    w = w3
    args = 'c2 c1'
  []

  [w33_res]
    type = ADSplitCHWRes
    variable = w3
    mob_name = L_33
  []

  [w32_res]
    type = ADSplitCHWRes
    variable = w3
    w = w2
    mob_name = L_32
  []

  [time_c2]
    type = ADCoupledTimeDerivative
    variable = w2
    v = c2
  []

  [time_c3]
    type = ADCoupledTimeDerivative
    variable = w3
    v = c3
  []
[]

[Materials]
# We need to change the length scale to units of nanometers.
# and change the energy scale to units of electron volts.
# The conversion from meters to nanometers is 1e+09
# The conversion from joules to electron volts is 6.24150934e+18

# Koyama 2004, Table 1: kappa_c = 1.0e-14 (J*m^2/mol)

  [gradient_coef_kappa_c]
    type = GenericFunctionMaterial
    prop_names = 'kappa_c'
    prop_values = '1.0e-14*6.24150934e+18*1e+09^2*1e-27' # eV*nm^2/mol
  []

  [constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' p    T    nm_m   eV_J            d      R                 Q_1     Q_2     Q_3     D01     D02     D03     bohr'
    prop_values = '0.4  913  1e+09  6.24150934e+18  1e-27  8.31446261815324  294000  308000  294000  1.0e-4  2.0e-5  1.0e-4  5.7883818012e-5'
  [] # bohr is Bohr magneton = eV/T

  [RT_clnc] # Koyama 2004, Second term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = RT_clnc
    material_property_names = 'T eV_J d R'
    args = 'c2 c3 c1'
    function = '1*eV_J*d*(R*T*(c1*log(c1)+c2*log(c2)+c3*log(c3)))'
    derivative_order = 2
  []

  [heat_of_mixing] # Koyama 2004, Third term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = E_G
    material_property_names = 'T eV_J d Lalpha_12(T) Lalpha_13(T) Lalpha_23(c2,c3,T)'
    function = 'eV_J*d*(Lalpha_12*c1*c2 + Lalpha_13*c1*c3 + Lalpha_23*c2*c3)'
    args = 'c2 c3 c1'
    derivative_order = 2
  []

  [magnetic_contribution_to_Gibbs] # Koyama 2004, Equation 2: ^mgG^alpha
    type = DerivativeParsedMaterial
    f_name = mg_G
    material_property_names = 'g(c1,c2,c3) beta(c1,c2,c3) T eV_J d R'
    function = '1*eV_J*d*(R*T*log(beta+1)*g)' #g is a tau function. Defined below
    args = 'c2 c3 c1'
    derivative_order = 2
  []

  [tau] # Koyama 2004, Equation 2: tau = T/(T^alpha_C)
    type = ParsedMaterial
    f_name = tau
    function = 'T/Tc'
    material_property_names = 'Tc(c1,c2,c3) T'
    args = 'c1 c2 c3'
    outputs = exodus
  []

  [tau_function] # Koyama 2004, Equation 2: f(tau)
  # Koyama does not explicitly define f(tau)
  # Expression for f(tau) was obtained by Xiong 2012, Equation 5
    type = ParsedMaterial
    f_name = g
    function='if(tau<1, 1-1/A*(79*tau^-1/(140*p)+474/497*(1/p-1)*(tau^3/6+tau^9/135+tau^15/600)), -1/A*(1/10*tau^-5+1/315*tau^-15+1/1500*tau^-25))'
    args = 'c1 c2 c3'
    material_property_names = 'p A tau(c1,c2,c3)'
  []

  [A] #Xiong 2012, Equation 6
  # A is a variable used to calculate tau function
    type = ParsedMaterial
    f_name = A
    function = '518/1125 + 11692/15975*(1/p-1)'
    material_property_names = 'p'
  []

  [Interaction_parameter_1_2] # Koyama 2004, Pg.2, 1st equation after Eq. 2: L^alpha_1,2
    type = DerivativeParsedMaterial
    f_name = Lalpha_12
    material_property_names = 'T'
    function = '20500 - 9.68*T'
    derivative_order = 2
  []

  [Interaction_parameter_1_3] # Koyama 2004, Pg.2, 2nd equation after Eq. 2: L^alpha_1,3
    type = DerivativeParsedMaterial
    f_name = Lalpha_13
    material_property_names = 'T'
    function = '-23669 + 103.9627*T - 12.7886*T*log(T)'
    derivative_order = 2
  []

  [Interaction_parameter_2_3] # Koyama 2004, Pg.2, 3rd equation after Eq. 2: L^alpha_2,3
    type = DerivativeParsedMaterial
    f_name = Lalpha_23
    material_property_names = 'T'
    function = '(24357 - 19.797*T) - 2010*(c3 - c2)'
    args = 'c2 c3'
    derivative_order = 2
  []

  [Curie_Temperature] # Koyama 2004, Pg.2, 4th equation after Eq. 2: T^alpha_C
    type = ParsedMaterial
    f_name = Tc
    function = '1043*c1 - 311.5*c2 + 1450*c3 + (1650 + 550*(c2-c1))*c1*c2 + 590*c1*c3'
    args = 'c1 c2 c3'
    outputs = exodus
  []

  [atomic_magnetic_moment] # Koyama 2004, Pg.2, 5th equation after Eq. 2: beta^alpha
    type = ParsedMaterial
    f_name = beta
    function = '(2.22*c1 - 0.01*c2 + 1.35*c3 - 0.85*c1*c2 + (2.4127 + 0.2418*(c3-c1))*c1*c3)'
    args = 'c1 c2 c3'
  []

  [Intensity_of_magnetization] # Koyama 2004, 2nd Page (Between Eq 6 & 7)
    type = ParsedMaterial
    f_name = I
    function='if(tau>0.9, 2^(-(2+10*(tau-1))), 1-1/7*(5*tau^4+2*tau^20))'
    args = 'c1 c2 c3'
    material_property_names = 'tau(c1,c2,c3)'
  []

  [magnetization_moment] # Koyama 2004, Equation 7
    type = ParsedMaterial
    f_name = Ms
    function = 'bohr*beta*I'
    args = 'c1 c2 c3'
    material_property_names = 'bohr beta(c1,c2,c3) I(c1,c2,c3)'
    outputs = exodus
  []

  [self_diffusion_coef_1] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_1
    material_property_names = 'D01 Q_1 R T'
    function = 'D01*exp(-Q_1/(R*T))'
    derivative_order = 1
  []

  [self_diffusion_coef_2] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_2
    material_property_names = 'D02 Q_2 R T'
    function = 'D02*exp(-Q_2/(R*T))'
    derivative_order = 1
  []

  [self_diffusion_coef_3] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_3
    material_property_names = 'D03 Q_3 R T'
    function = 'D03*exp(-Q_3/(R*T))'
    derivative_order = 1
  []

# Onsager coefficients: "Mobility_Lij":
  [Mobility_L22] # Koyama 2004, After Equation 8: L22
    type = DerivativeParsedMaterial
    f_name = L_22
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*c2*D_1 + (1-c2)^2*D_2 + c2*c3*D_3)*c2/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [Mobility_L23] # Koyama 2004, After Equation 8: L23
    type = DerivativeParsedMaterial
    f_name = L_23
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*D_1 - (1-c2)*D_2 - (1-c3)*D_3)*c2*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [Mobility_L32] # Koyama 2004, After Equation 8: L32 = L23
    type = DerivativeParsedMaterial
    f_name = L_32
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*D_1 - (1-c2)*D_2 - (1-c3)*D_3)*c2*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [Mobility_L33] # Koyama 2004, After Equation 8: L33
    type = DerivativeParsedMaterial
    f_name = L_33
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*c3*D_1 + c2*c3*D_2 + (1-c3)^2*D_3)*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [G_system] # Total free energy of the system
    type = DerivativeSumMaterial
    f_name = F_total
    sum_materials = 'E_G RT_clnc mg_G'
    args = 'c1 c2 c3'
    outputs = exodus
    derivative_order = 2
  []
[]

[Preconditioning]
  [coupled]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  steady_state_detection = true
  solve_type = PJFNK
  automatic_scaling = true
  l_max_its = 30
  l_tol = 1e-6
  nl_max_its = 50
  nl_abs_tol = 1e-8
  end_time = 300   # 1h
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type
                         -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly
                         ilu          1'

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    cutback_factor = 0.8
    growth_factor = 1.2
    optimal_iterations = 9
  []

  [Adaptivity]
    coarsen_fraction = 0.1
    refine_fraction = 0.7
    max_h_level = 2
  []
[]

[Outputs]
  exodus = true
[]
