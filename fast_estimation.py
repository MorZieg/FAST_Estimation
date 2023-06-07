# FAST Estimation v1.1 - GPLv3
# Moritz O. Ziegler, mziegler@gfz-potsdam.de
# Manual:   https://doi.org/10.48440/wsm.2023.001
# Download: https://github.com/MorZieg/FAST_Estimation
###############################################################################################
#
# Python Fast Automatic Stress Tensor Estimation v1.1 is a Python 3 tool of the FAST suite.
# It supports Moose and Abaqus solver, PyTecplot, GeoStress and GeoStressCmd
# (https://doi.org/10.5880/wsm.2020.001), and runs on Windows and Linux
# systems. It can be run as stand-alone script or called from another script.
#
# The basic principle is explained in the FAST Estimation Manual referenced above. The
# additionally required information is provided in this file. For simple usage refer to
# the user input section only. Additional comments are provided throughout the script.
#
# Important information:
# When using Tecplot Macros (pytecplot = 'off') the three steps with the different boundary
# conditions are expected to be loaded in Tecplot. Make sure that the folder variable is set
# correctly and that the subdirectory 'data' exists in the working directory.
# 
# When using PyTecplot:
#   Abaqus:  The native variables stress variables (XX Stress etc.) or SHmax & Shmin are
#        expected.
#      On Windows a *.odb file is expected. (May have to be loaded before to be
#         converted to a readable format.
#      On Linux a *.fil file is expected.
#   Moose:   Three *.dat Tecplot files provided as output by Moose are required. After
#        the name they are expected to be named *_0001.dat, *_0002.dat, and *_0003.dat.
#      The native stress variables (stress_xx etc.) or SHmax & Shmin are expected.
###############################################################################################

#  User input:
solver = 'abaqus' # 'moose'
pytecplot = 'on' # 'off'

#  The stress components that should be estimated
#stress_vars = ['XX Stress','YY Stress','ZZ Stress','XY Stress','YZ Stress','ZX Stress']
stress_vars = ['SHmax','Shmin','Sv']

#  The name of the test file with three boundary condition scenarios.
name = 'test_scenarios'

#  The current folder.
folder = 'C:\\Users\\Documents\\Software\\FAST_Estimation'

#  The test boundary conditions.
bcs = [[4, -4],[2, -5],[4, -3]]

# The new boundary conditions that should be evaluated
#bc_eval = 'bc_estimation.csv'
bc_eval = [[ -4.4, 2],[-6, 4.1]]

#  Locations (X Y Z) to be evaluated
#loc = 'locations.csv'
loc = [[3500, 3500, -2900],[8000, 7000, -900]]

###############################################################################################
def main(loc,bc_eval,stress_vars,bcs,name,solver,pytecplot):
  import os.path
  import numpy as np
  import tecplot

  print('Running FAST Estimation v1.1')
  
  # Load locations if specified in file.
  if type(loc) == str:
    loc = load_loc(loc)

  # Load boundary conditions if specified in file.
  if type(bc_eval) == str:
    bc_eval = load_bc(bc_eval)

  # Load reference data from three test scenarios
  if pytecplot == 'on':
    print('Running PyTecplot Version '+tecplot.__version__)
    # Load files if PyTecplot is used
    if not os.path.exists((name+'.plt')):
      if solver == 'abaqus':
        load_abq(name)
      elif solver == 'moose':
        load_mse(name)
    else:
      print('Using existing *.plt file')
      
    moss = extract_tp(name,solver,loc,stress_vars)

  elif pytecplot == 'off':
    # If PyTecplot is NOT used, create a macro
    write_macro(loc,stress_vars,name,folder)
    
    # A break appears here in order to execute the macro in Tecplot.
    print('Execute Macro in Tecplot...')
    input("...then press Enter to continue...")

    # Load variables from .csv file.
    moss = []
    for i in range(len(stress_vars)):
      moss.append(load_csv(name,len(loc),stress_vars[i]))
        
  else:
    print('ERROR! Indicate if you are using Pytecplot or not.')

  # The actual estimation starts here!
  es = []
  for j in range(len(bc_eval)):
    es.append([bc_eval[j]])
    shazi = 0
    for i in range(len(loc)):
      temp = []
      for k in range(len(stress_vars)):        
        temp.append(solve(moss[k][i],bcs,bc_eval[j][0],bc_eval[j][1]))

      es[j].append(temp)

    es[j][0].append(stress_rotation(es[j],stress_vars))
  
  return es

###############################################################################################
###############################################################################################
def write_macro(coords,stress_vars,name,folder):
  # A Tecplot macro to export the stress state at the desired locations is written.
  import numpy as np
  import datetime

  now = datetime.datetime.now()
  timestamp = str(now.strftime("%Y%m%d_%H%M%S"))
  fid =  open(('macro_FAST_Estimation_'+timestamp+'.mcr'),'w')
  ########## WARNING! ############
  # The following line contains the Tecplot Macro header.
  # The number may need to be adapted for different versions of Tecplot.
  fid.write('#!MC 1410\n\n')
  ################################
  output = ''
  for i in range(len(stress_vars)):
    # Find the number Tecplot assigned to the variables according to their names.
    fid.write('$!GETVARNUMBYNAME |VAR%i|\nNAME = "%s"\n' % (i, stress_vars[i]))
    output = output + '|VAR' + str(i) + '|,'
  output = output[0:-1]

  # Create 1D zones at the locations. At each point 3 zones are created.
  for i in range(len(coords)):
    for j in range(3):
      fid.write('$!CREATERECTANGULARZONE\nIMAX = 1\nJMAX = 1\nKMAX = 1\n')
      fid.write('X1 = %i\nY1 = %i\nZ1 = %i\n' % (coords[i][0],coords[i][1],coords[i][2]))
      fid.write('X2 = %i\nY2 = %i\nZ2 = %i\n\n' % (coords[i][0],coords[i][1],coords[i][2]))
  
  # Interpolate stress state to newly created 1D zones. From each boundary condition
  # scenario at each location the stress state is interpolated to a 1D zone.
  for i in range(len(coords)):
    for x,j in enumerate([1,2,3]):
      fid.write('$!LINEARINTERPOLATE\nSOURCEZONES =  [%i]\n' % j)
      fid.write('DESTINATIONZONE = %i\nVARLIST =  [%s]\nLINEARINTERPCONST = 0\nLINEARINTERPMODE = DONTCHANGE\n\n' % ((3+(i*3)+j),output))
  
  # Export the three stress components to individually named files in the specified folder.
  for i in range(len(stress_vars)):
    inst = 'VAR' + str(i)
    id = '_' + stress_vars[i]
    
    fid.write('$!EXTENDEDCOMMAND\nCOMMANDPROCESSORID = \'excsv\'\n')
    fid.write('COMMAND = \'FrOp=1:ZnCount=%i:ZnList=[%i-%i]:VarCount=1:VarList=[|%s|]:ValSep=",":FNAME="%s\\data\\%s.csv"\'\n\n' % (len(coords)*3,4,len(coords)*3+3,inst,folder,(name+id)))
  
  # Delet the 1D zones.
  fid.write('$!DELETEZONES [%i-%i]' % (4,len(coords)*3+3))
  
  fid.close()

###############################################################################################
def load_abq(name):
  # Only required if PyTecplot is used.
  # Load Abaqus output file and save it as *.plt file. Different input files
  # (and syntax) required depending on operating system.
  import tecplot as tp
  import platform

  print('Loading *.odb file')
  if platform.system() == 'Linux':
    tp.macro.execute_command("""$!ReadDataSet  '\"StandardSyntax\" \"1.0\" \"FEALoaderVersion\" \"446\" \"FILENAME_File\" \"%s.fil\" \"AutoAssignStrandIDs\" \"Yes\"'\n"""
    """  DataSetReader = 'ABAQUS .fil Data (FEA)'""" % name)
    
  elif platform.system() == 'Windows':
    tp.macro.execute_command("""$!ReadDataSet  '\"StandardSyntax\" \"1.0\" \"FEALoaderVersion\" \"446\" \"FILENAME_File\" \"%s.odb\" \"AutoAssignStrandIDs\" \"Yes\"'"""
    """DataSetReader = 'ABAQUS Output Database (FEA)'""" % name)
    
  else:
    print ('Platform/OS not recognized')
  
  # Save as *.plt file.
  tp.data.save_tecplot_plt('%s.plt' % name, include_text=False, include_geom=False, include_data_share_linkage=True)
  
###############################################################################################
def load_mse(name):
  # Only required in PyTecplot is used.
  # Load Moose output files consecutively and save as *.plt file.
  import tecplot as tp

  print('Loading Moose output file')
  # Load first solver output file.
  tp.macro.execute_command("""$!ReadDataSet  '\"%s_0001.dat\" '
    ReadDataOption = New
    ResetStyle = No
    VarLoadMode = ByName
    AssignStrandIDs = Yes
    VarNameList = '\"x\" \"y\" \"z\" \"stress_xx\" \"stress_yy\" \"stress_zz\" \"stress_xy\" \"stress_yz\" \"stress_zx\" \"disp_x\" \"disp_y\" \"disp_z\"'""" % name)
  
  # Append second solver output files.
  tp.macro.execute_command("""$!ReadDataSet  '\"%s_0002.dat\" '
    ReadDataOption = Append
    ResetStyle = No
    VarLoadMode = ByName
    AssignStrandIDs = Yes
    VarNameList = '\"x\" \"y\" \"z\" \"stress_xx\" \"stress_yy\" \"stress_zz\" \"stress_xy\" \"stress_yz\" \"stress_zx\" \"disp_x\" \"disp_y\" \"disp_z\"'""" % name)

  # Append third solver output files.
  tp.macro.execute_command("""$!ReadDataSet  '\"%s_0003.dat\" '
    ReadDataOption = Append
    ResetStyle = No
    VarLoadMode = ByName
    AssignStrandIDs = Yes
    VarNameList = '\"x\" \"y\" \"z\" \"stress_xx\" \"stress_yy\" \"stress_zz\" \"stress_xy\" \"stress_yz\" \"stress_zx\" \"disp_x\" \"disp_y\" \"disp_z\"'""" % name)
  
  # Save as *.plt file.
  tp.data.save_tecplot_plt('%s.plt' % name, include_text=False, include_geom=False, include_data_share_linkage=True)
  
###############################################################################################
def extract_tp(name,solver,loc,stress_vars):
  # Read the *.plt file created by load_abq or load_mse. If the stress variables
  # SHmax and Shmin are desired, GeoStressCmd is run to compute them. Then the
  # modelled stress state (moss) is extracted at the according locations using strextract.
  import tecplot as tp
  import numpy as np
  import platform
  
  model = tp.data.load_tecplot('%s.plt' % name)
  
  if platform.system() == 'Linux' and solver == 'abaqus':
    # Convert the cell centered stress tensor variables from the *.fil file to nodal variables.
    cell2nodal(model)
  
  if solver == 'moose':
    # Rename variables to those that can be read by GeoStress
    rnm_vrbls()
     
  # Run GeoStressCmd if a reduced stress tensor is desired.
  if stress_vars[0] == 'SHmax' or stress_vars[0] == 'Shmin':
    print ('Running GeoStressCmd...')
    CommandString = 'Stress, 90.0, 0.0, 0.0, 0.0, -1.0E-6, 0, 0, 0, 0, 0, 1, 0, 0, 0'
    tp.macro.execute_extended_command("GeoStressCmd",CommandString)
    print ('Sucessfull!')
  
  # Extract the variables at the acording locations.
  moss = []
  for j in range(len(stress_vars)):
    temp = np.zeros((len(loc),3))
    for i in range(len(loc)):
      for k in range(3):
        temp[i][k] = strextract(loc[i][0],loc[i][1],loc[i][2],k,model,stress_vars[j])

    moss.append(temp)
  
  return moss

###############################################################################################
def solve(moss,bcs,bcnx,bcny):
  import numpy as np
  # moss: Data from three test scenarios for the stress component that is solved.
  # bcs: Displacement in x’ and y' direction prescribed at different test scenarios.
  # bcnx: Boundary conditions for estimated stress state in x' direction.
  # bcny: Boundary conditions for estimated stress state in y' direction.

  # Setup parameter form of plain
  p1 = np.array([bcs[0][0], bcs[0][1], moss[0]])
  p2 = np.array([bcs[1][0], bcs[1][1], moss[1]])
  p3 = np.array([bcs[2][0], bcs[2][1], moss[2]])

  p = p1
  u = p1 - p2
  v = p1 - p3

  # Check if the data is linearly independent
  independence_check(u,v)
  
  # normal form / vector
  n = np.cross(u,v)

  # coordinate form
  d = p[0] * n[0] + p[1] * n[1] + p[2] * n[2]

  stress = ( d - n[0] * bcnx - n[1] * bcny ) / n[2]

  return stress

###############################################################################################
###############################################################################################
def cell2nodal(model):
  import tecplot as tp
  from tecplot.constant import ValueLocation
  
  # Create new stress tensor variables at nodes.
  tp.data.operate.execute_equation('{XX Stress node} = {XX Stress}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{YY Stress node} = {YY Stress}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{ZZ Stress node} = {ZZ Stress}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{XY Stress node} = {XY Stress}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{YZ Stress node} = {YZ Stress}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{ZX Stress node} = {ZX Stress}', value_location=ValueLocation.Nodal)
  
  # Delete the original (cell-centered) stress tensor variables.
  model.delete_variables(model.variable('XX Stress'),model.variable('YY Stress'),model.variable('ZZ Stress'))
  model.delete_variables(model.variable('XY Stress'),model.variable('YZ Stress'),model.variable('ZX Stress'))
  
  # Create stress tensor variables at nodes.
  tp.data.operate.execute_equation('{XX Stress} = {XX Stress node}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{YY Stress} = {YY Stress node}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{ZZ Stress} = {ZZ Stress node}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{XY Stress} = {XY Stress node}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{YZ Stress} = {YZ Stress node}', value_location=ValueLocation.Nodal)
  tp.data.operate.execute_equation('{ZX Stress} = {ZX Stress node}', value_location=ValueLocation.Nodal)
  
  # Delete the temporal stress tensor variables.
  model.delete_variables(model.variable('XX Stress node'),model.variable('YY Stress node'),model.variable('ZZ Stress node'))
  model.delete_variables(model.variable('XY Stress node'),model.variable('YZ Stress node'),model.variable('ZX Stress node'))
  
###############################################################################################
def independence_check(u,v):
  # Check whether the test boundary scenarios are correctly chosen and independent.
  for i in range(len(v)):
    if v[i] == 0:
      v[i] = 0.0001
  
  test = u/v;

  if test[1] == test[2] == test[3]:
    print('ERROR! Planes are linearly dependent. Choose different initial boundary conditions')

###############################################################################################
def load_csv(name,leloc,stress_vars):
  # Load the files/stress components that were exported from Tecplot with the macro. 
  import csv
  import numpy as np
  
  # Read the files that contain the variables.
  temp = []
  moss = []
  with open('data/' + name + '_' + stress_vars + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      if row != []:
        temp.append(float(row[0]))
       
  # Sort the variables according to locations.
  for i in range(leloc):
    temp2 = []
    for j in range(3):
      temp2.append(temp[j+i*3])
    moss.append(temp2)
    
  # The variables are converted to a numpy array.
  moss = np.array(moss)
  
  return moss

###############################################################################################
def load_bc(bc):
  # Load the boundary conditions used to estimate the stress state from a .csv file.
  import csv
  import numpy as np
  
  bc_eval = []
  with open(bc) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      if row != []:
        bc_eval.append([float(row[0]),float(row[1])])
        
  return bc_eval

###############################################################################################
def load_loc(loc):
  # Load the locations where the stress state is estimated from a .csv file.
  import csv
  import numpy as np
  
  coord = []
  with open(loc) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
      if row != []:
        coord.append([float(row[0]),float(row[1]),float(row[2])])
        
  return coord

###############################################################################################
def rnm_vrbls():
  import tecplot as tp
  
  # Assign solution time.
  tp.macro.execute_extended_command(command_processor_id='Strand Editor',command='ZoneSet=1-3;AssignStrands=TRUE;StrandValue=1;AssignSolutionTime=TRUE;TimeValue=1;DeltaValue=1;TimeOption=ConstantDelta;')
  
  # Generate variables according to the naming convetion of GeoStress.
  tp.data.operate.execute_equation(equation='{XX Stress} = {stress_xx}')
  tp.data.operate.execute_equation(equation='{YY Stress} = {stress_yy}')
  tp.data.operate.execute_equation(equation='{ZZ Stress} = {stress_zz}')
  tp.data.operate.execute_equation(equation='{XY Stress} = {stress_xy}')
  tp.data.operate.execute_equation(equation='{YZ Stress} = {stress_yz}')
  tp.data.operate.execute_equation(equation='{ZX Stress} = {stress_zx}')
  
  tp.data.operate.execute_equation(equation='{X Displacement} = {disp_x}')
  tp.data.operate.execute_equation(equation='{Y Displacement} = {disp_y}')
  tp.data.operate.execute_equation(equation='{Z Displacement} = {disp_z}')

  tp.data.save_tecplot_plt('%s.plt' % 'test_moose', include_text=False, include_geom=False, include_data_share_linkage=True)

###############################################################################################
def strextract(x,y,z,zone,model,comp):
  # Function that extracts specified variables at a certain location and zone.
  import tecplot as tp
  
  tecvar = model.variable(comp).index
  result = tp.data.query.probe_at_position(x,y,z,zones=[zone])
  stress = result[0][tecvar]
  
  return stress

###############################################################################################
def stress_rotation(temp,stress_vars):
  # Check whether a stress rotation takes place, i.e. if Shmin > SHmax.
  shazi = 0

  if stress_vars[0] == 'SHmax' and stress_vars[1] == 'Shmin':
    for i in range(1,len(temp)):
      if temp[i][0] > temp[i][1]:
        shazi = shazi + 1
        
    if shazi == 0:
      out = 'SHmax azimuth remains.'
    elif shazi > 0:
      out = 'ATTENTION! SHmax azimuth changes at ' + str(shazi) + ' locations.'

  else:
    out = 'No stress rotation detection was carried out.'
    

  return out

###############################################################################################
def write_output(es,stress_vars):
  # Write the estimated stress states to a *.dat file
  import datetime

  now = datetime.datetime.now()
  timestamp = str(now.strftime("%Y%m%d_%H%M%S"))
  fid =  open(('estimated_stress_states_'+timestamp+'.dat'),'w')
  fid.write('Stress states estimated by FAST Estimation v1.0\n')
  fid.write('Stress components: ')
  for i in range(len(stress_vars)):
    fid.write('%s ' % stress_vars[i])
  fid.write('\n')
  
  for i in range(len(es)):
    fid.write('# Boundary Condition Scenario # %i\n' % (i+1) )
    fid.write('%f, %f' % (es[i][0][0], es[i][0][1]) )
    if len(es[i][0]) > 2:
      fid.write(', %s' % es[i][0][2])
    fid.write('\n# Estimated Stress States\n')
    for j in range(1,len(es[i])):
      for k in range(len(es[0][1])):
        fid.write('%.4f, ' % es[i][j][k] )
      fid.write('\n')
    fid.write('#\n')
  fid.close()

###############################################################################################
###############################################################################################
if __name__ == '__main__':
  es = main(loc,bc_eval,stress_vars,bcs,name,solver,pytecplot)
  
  # If the function is run as a stand-alone script the number of estimated stress states are
  # printed to the screen and a results file is written.
  print(str(len(es)) + ' stress scenarios estimated at')
  print(str(len(es[0])-1) + ' individual locations.')
  write_output(es,stress_vars)
  print('Results written to file.')
###############################################################################################
