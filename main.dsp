# Microsoft Developer Studio Project File - Name="main" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=main - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "main.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "main.mak" CFG="main - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "main - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "main - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "main - Win32 Release"
# Name "main - Win32 Debug"
# Begin Source File

SOURCE=.\block.f90
DEP_F90_BLOCK=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\run_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\blockgmres.f90
DEP_F90_BLOCKG=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\run_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\correct.f90
# End Source File
# Begin Source File

SOURCE=.\deformation_gradient.f90
DEP_F90_DEFOR=\
	".\fluid_variables.mod"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\delta_nonuniform.f90
DEP_F90_DELTA=\
	".\fluid_variables.mod"\
	".\solid_variables.mod"\
	".\vol2d3n.fi"\
	".\vol2d4n.fi"\
	".\vol3d4n.fi"\
	
# End Source File
# Begin Source File

SOURCE=.\echoinput.f90
DEP_F90_ECHOI=\
	".\fluid_variables.mod"\
	".\run_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ensight_output.f90
DEP_F90_ENSIG=\
	".\fluid_variables.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\equal.f90
# End Source File
# Begin Source File

SOURCE=.\error.f90
# End Source File
# Begin Source File

SOURCE=.\facemap.f90
DEP_F90_FACEM=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\fclear.f90
# End Source File
# Begin Source File

SOURCE=.\fluid_variables.f90
# End Source File
# Begin Source File

SOURCE=.\form.f90
DEP_F90_FORM_=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\formm.f90
DEP_F90_FORMM=\
	".\fluid_variables.mod"\
	".\run_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\gaussj.f90
# End Source File
# Begin Source File

SOURCE=.\gjinv.f90
# End Source File
# Begin Source File

SOURCE=.\global_constants.f90
# End Source File
# Begin Source File

SOURCE=.\global_simulation_parameter.f90
# End Source File
# Begin Source File

SOURCE=.\gmres.f90
DEP_F90_GMRES=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\gmres_fluid_ale.f90
DEP_F90_GMRES_=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\run_variables.mod"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\gmres_fluid_mesh.f90
DEP_F90_GMRES_F=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\hg.f90
DEP_F90_HG_F9=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\hydro.f90
DEP_F90_HYDRO=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\hypo.f90
DEP_F90_HYPO_=\
	".\delta_nonuniform.mod"\
	".\ensight_output.mod"\
	".\fluid_variables.mod"\
	".\form.mod"\
	".\global_constants.mod"\
	".\global_simulation_parameter.mod"\
	".\hypo_declaration_fluid.fi"\
	".\hypo_declaration_solid.fi"\
	".\hypo_fluid_solver.fi"\
	".\hypo_prepare_fluid.fi"\
	".\hypo_prepare_solid.fi"\
	".\hypo_restart_file_check.fi"\
	".\hypo_restart_read.fi"\
	".\hypo_restart_write.fi"\
	".\hypo_write_output.fi"\
	".\meshgen_fluid.mod"\
	".\meshgen_solid.mod"\
	".\r_common.mod"\
	".\run_variables.mod"\
	".\solid_fem_BC.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\hypo_declaration_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_declaration_solid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_fluid_solver.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_fluid_solver_mesh.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_prepare_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_prepare_solid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_file_check.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_read.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_write.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_write_output.fi
# End Source File
# Begin Source File

SOURCE=.\initialize.f90
DEP_F90_INITI=\
	".\fluid_variables.mod"\
	".\run_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\lenght.f90
DEP_F90_LENGH=\
	".\fluid_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\length.f90
DEP_F90_LENGT=\
	".\fluid_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\parseinput.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\meshgen_fluid.f90
DEP_F90_MESHG=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\meshgen_solid.f90
DEP_F90_MESHGE=\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nondimension.f90
DEP_F90_NONDI=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\norm.f90
# End Source File
# Begin Source File

SOURCE=.\normal.f90
DEP_F90_NORMA=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\parseinput.f90
DEP_F90_PARSE=\
	".\delta_nonuniform.mod"\
	".\fluid_variables.mod"\
	".\global_simulation_parameter.mod"\
	".\meshgen_solid.mod"\
	".\r_common.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\quad2d3n.f90
# End Source File
# Begin Source File

SOURCE=.\quad2d4n.f90
# End Source File
# Begin Source File

SOURCE=.\quad3d4n.f90
# End Source File
# Begin Source File

SOURCE=.\quad3d8n.f90
# End Source File
# Begin Source File

SOURCE=.\r_bdpd_curr.f90
DEP_F90_R_BDP=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_bdpd_init.f90
DEP_F90_R_BDPD=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_common.f90
# End Source File
# Begin Source File

SOURCE=.\r_element.f90
DEP_F90_R_ELE=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_jacob.f90
DEP_F90_R_JAC=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_load.f90
DEP_F90_R_LOA=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_nodalf.f90
DEP_F90_R_NOD=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sboc.f90
DEP_F90_R_SBO=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sbpress.f90
DEP_F90_R_SBP=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scal.f90
DEP_F90_R_SCA=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scauchy.f90
DEP_F90_R_SCAU=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_smaterj.f90
DEP_F90_R_SMA=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spiola.f90
DEP_F90_R_SPI=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spiola_elastic.f90
DEP_F90_R_SPIO=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spiola_viscous.f90
DEP_F90_R_SPIOL=\
	".\fluid_variables.mod"\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spress.f90
DEP_F90_R_SPR=\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sreadinit.f90
DEP_F90_R_SRE=\
	".\r_common.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sstif.f90
DEP_F90_R_SST=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sstrain.f90
DEP_F90_R_SSTR=\
	".\r_common.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_stang.f90
DEP_F90_R_STA=\
	".\r_common.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_stoxc.f90
DEP_F90_R_STO=\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_timefun.f90
DEP_F90_R_TIM=\
	".\global_constants.mod"\
	".\r_common.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\read.f90
# End Source File
# Begin Source File

SOURCE=.\rkpmshape2d.f90
# End Source File
# Begin Source File

SOURCE=.\rkpmshape3d.f90
# End Source File
# Begin Source File

SOURCE=.\run_variables.f90
# End Source File
# Begin Source File

SOURCE=.\set.f90
DEP_F90_SET_F=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\shape.f90
DEP_F90_SHAPE=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\sharp.f90
DEP_F90_SHARP=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\solid_fem_BC.f90
DEP_F90_SOLID=\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\solid_solver.f90
DEP_F90_SOLID_=\
	".\r_common.mod"\
	".\solid_fem_BC.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\solid_update.f90
DEP_F90_SOLID_U=\
	".\r_common.mod"\
	".\run_variables.mod"\
	".\solid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\solid_variables.f90
# End Source File
# Begin Source File

SOURCE=.\update.f90
DEP_F90_UPDAT=\
	".\fluid_variables.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\velocity.f90
DEP_F90_VELOC=\
	".\fluid_variables.mod"\
	".\global_constants.mod"\
	".\run_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\vol.f90
DEP_F90_VOL_F=\
	".\fluid_variables.mod"\
	".\sh2d3n.h"\
	".\sh2d4n.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\vol2d3n.fi
# End Source File
# Begin Source File

SOURCE=.\vol2d4n.fi
# End Source File
# Begin Source File

SOURCE=.\vol3d4n.fi
# End Source File
# End Target
# End Project
