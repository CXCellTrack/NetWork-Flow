@echo ��� ssvm train test �û���������
@echo off
rem ע�Ⲣ���û�������
rem ע��xxx=XXX�м�һ�����ܼӿո񣡷���ᱨ��

rem set regpath=HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment
set regpath=HKEY_CURRENT_USER\Environment

set evname=trainSegPath

set trainSegPath=E:\datasets\first_edition\training_datasets\N2DH-GOWT1\05_0-64_seg
reg add "%regpath%" /v %evname% /d %trainSegPath% /f

set evname=trainTrackPath

set trainTrackPath=E:\datasets\first_edition\training_datasets\N2DH-GOWT1\05_0-64_track

reg add "%regpath%" /v %evname% /d %trainTrackPath% /f



set evname=testSegPath

set testSegPath=E:\datasets\second_edition\training_datasets\Fluo-N2DH-SIM+\02_4-16_seg

reg add "%regpath%" /v %evname% /d %testSegPath% /f

set evname=testTrackPath

set testTrackPath=E:\datasets\second_edition\training_datasets\Fluo-N2DH-SIM+\02_4-16_track

reg add "%regpath%" /v %evname% /d %testTrackPath% /f






rem ɾ��ϵͳ���� reg delete "%regpath%" /v "%evname%"  /f

rem �������������������ˢ��ע���

rem taskkill /f /im explorer.exe && explorer.exe 
 
rem ��ͣ������ʾ���밴�����������

ping localhost -n 3 >nul

exit

rem pause>nul


