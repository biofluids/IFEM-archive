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
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

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
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "main - Win32 Release"
# Name "main - Win32 Debug"
# Begin Source File

SOURCE=.\add.f
# End Source File
# Begin Source File

SOURCE=.\block.f
# End Source File
# Begin Source File

SOURCE=.\blockgmres.f
# End Source File
# Begin Source File

SOURCE=.\blockgmresm.f
# End Source File
# Begin Source File

SOURCE=.\blockm.f
# End Source File
# Begin Source File

SOURCE=.\cut.f
# End Source File
# Begin Source File

SOURCE=.\dataout.c
# End Source File
# Begin Source File

SOURCE=.\defgrad.f
# End Source File
# Begin Source File

SOURCE=.\dim.f
# End Source File
# Begin Source File

SOURCE=.\disk.f
# End Source File
# Begin Source File

SOURCE=.\dyno.f
# End Source File
# Begin Source File

SOURCE=.\echoinput.f
# End Source File
# Begin Source File

SOURCE=.\equal.f
# End Source File
# Begin Source File

SOURCE=.\error.f
# End Source File
# Begin Source File

SOURCE=.\ewdio.c
# End Source File
# Begin Source File

SOURCE=.\ewdmem.c
# End Source File
# Begin Source File

SOURCE=.\fcc.f
# End Source File
# Begin Source File

SOURCE=.\fclear.f
# End Source File
# Begin Source File

SOURCE=.\form.f
# End Source File
# Begin Source File

SOURCE=.\formm.f
# End Source File
# Begin Source File

SOURCE=.\get.f
# End Source File
# Begin Source File

SOURCE=.\getint.c
# End Source File
# Begin Source File

SOURCE=.\getkey.c
# End Source File
# Begin Source File

SOURCE=.\getreal.c
# End Source File
# Begin Source File

SOURCE=.\getstr.c
# End Source File
# Begin Source File

SOURCE=.\getstrn.c
# End Source File
# Begin Source File

SOURCE=.\global.h
# End Source File
# Begin Source File

SOURCE=.\gmres.f
# End Source File
# Begin Source File

SOURCE=.\gmresm.f
# End Source File
# Begin Source File

SOURCE=.\hg.f
# End Source File
# Begin Source File

SOURCE=.\hydro.f
# End Source File
# Begin Source File

SOURCE=.\hypo.f
# End Source File
# Begin Source File

SOURCE=.\iei2cray.c
# End Source File
# Begin Source File

SOURCE=.\initialize.f
# End Source File
# Begin Source File

SOURCE=.\isatty.c
# End Source File
# Begin Source File

SOURCE=.\lenght.f
# End Source File
# Begin Source File

SOURCE=.\liftdrag.f
# End Source File
# Begin Source File

SOURCE=.\locate.f
# End Source File
# Begin Source File

SOURCE=.\lump.f
# End Source File
# Begin Source File

SOURCE=.\main.f
DEP_F90_MAIN_=\
	".\global.h"\
	
NODEP_F90_MAIN_=\
	".\Debug\mpif.h"\
	
# End Source File
# Begin Source File

SOURCE=.\malloc.h
# End Source File
# Begin Source File

SOURCE=.\mesh.info
# End Source File
# Begin Source File

SOURCE=.\meshgen.f
# End Source File
# Begin Source File

SOURCE=.\new.h
# End Source File
# Begin Source File

SOURCE=.\norm.f
# End Source File
# Begin Source File

SOURCE=.\parseinput.f
# End Source File
# Begin Source File

SOURCE=.\quad2d3n.f
# End Source File
# Begin Source File

SOURCE=.\quad2d4n.f
# End Source File
# Begin Source File

SOURCE=.\quad3d4n.f
# End Source File
# Begin Source File

SOURCE=.\quad3d8n.f
# End Source File
# Begin Source File

SOURCE=.\readdata.c
# End Source File
# Begin Source File

SOURCE=.\set.f
# End Source File
# Begin Source File

SOURCE=.\sh3d4n.h
# End Source File
# Begin Source File

SOURCE=.\sh3d8n.h
# End Source File
# Begin Source File

SOURCE=.\shape.f
# End Source File
# Begin Source File

SOURCE=.\sharp.f
# End Source File
# Begin Source File

SOURCE=.\solerror.f
# End Source File
# Begin Source File

SOURCE=.\update.f
# End Source File
# Begin Source File

SOURCE=.\updatenode.f
# End Source File
# Begin Source File

SOURCE=.\updatex.f
# End Source File
# Begin Source File

SOURCE=.\velocity.f
# End Source File
# Begin Source File

SOURCE=.\vol.f
# End Source File
# End Target
# End Project
