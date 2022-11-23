    
--[[
Sequence being followed

A. CONFIGURATION
1. Configuring Master from SOP till RF Init
2. Configuring Slave (i) sequentially from SOP till RF Init . i = 1, 2, 3 
3. Configuring Master from Data Path Config till Datapath Configuration.
4. Configuring Slaves from Data Path Config till Datapath Configuration.

NOTE:
Update the following in the script accordingly before running
1. metaImage F/W path on line 32
2. TDA Host Board IP Address on line 39
--]]

--------------------------------User Constants---------------------------------------
    
dev_list                 =    {1, 2, 4, 8}      -- Device map
RadarDevice              =    {1, 1, 1, 1}      -- {dev1, dev2, dev3, dev4}, 1: Enable, 0: Disable
                                                
cascade_mode_list        =    {1, 2, 2, 2}      -- 0: Single chip, 1: Master, 2: Slave

-- F/W Download Path

-- Uncomment the next line if you wish to pop-up a dialog box to select the firmware image file
-- Otherwise, hardcode the path to the firmware metaimage below
-- By default, the firmware filename is: xwr22xx_metaImage.bin
metaImagePath            =   "C:\\ti\\mmwave_dfp_02_02_03_01\\firmware\\xwr22xx_metaImage.bin"

-- IP Address for the TDA2 Host Board
-- Change this accordingly for your setup

TDA_IPAddress            =   "192.168.33.180"

-- Device map of all the devices to be enabled by TDA
-- 1 - master ; 2- slave1 ; 4 - slave2 ; 8 - slave3

deviceMapTDA             =   15                 -- all devices enabled

-------------------------------  API Configuration -----------------------------------

-- 1. Connection to TDA. 2. Selecting Cascade/Single Chip.  3. Selecting 2-chip/4-chip


WriteToLog("Setting up Studio for Cascade started..\n", "blue")

if (0 == ar1.ConnectTDA(TDA_IPAddress, 5001, deviceMapTDA)) then
    WriteToLog("ConnectTDA Successful\n", "green")
else
    WriteToLog("ConnectTDA Failed\n", "red")
    return -1
end

if (0 == ar1.selectCascadeMode(1)) then
    WriteToLog("selectCascadeMode Successful\n", "green")
else
    WriteToLog("selectCascadeMode Failed\n", "red")
    return -1
end

WriteToLog("Setting up Studio for Cascade ended..\n", "blue")
      


for i=1,table.getn(RadarDevice) do 

    local status    =    0        
    ---------------------------Configuration of AWR Devices------------------------
    
    if ((RadarDevice[1]==1) and (RadarDevice[i]==1)) then
    
        WriteToLog("Device"..i.." Configuration Started... \n", "blue")
        
        ---------------------------Connection Tab-------------------------
        
        -- SOP Mode Configuration
        if (0 == ar1.SOPControl_mult(dev_list[i], 4)) then
            WriteToLog("Device "..i.." : SOP Reset Successful\n", "green")
        else
            WriteToLog("Device "..i.." : SOP Reset Failed\n", "red")
            return -1
        end
                
        -- SPI Connect
        if (i==1) then
            status    =    ar1.PowerOn_mult(dev_list[i], 0, 1000, 0, 0)
        else
            status    =    ar1.AddDevice(dev_list[i])
        end        
        if (0 == status) then
            WriteToLog("Device "..i.." : SPI Connection Successful\n", "green")
        else
            WriteToLog("Device "..i.." : SPI Connection Failed\n", "red")
            return -1
        end
        
        -- Firmware Download. (SOP 4 - MetaImage)
        if (0 == ar1.DownloadBssFwOvSPI_mult(dev_list[i], metaImagePath)) then
            WriteToLog("Device "..i.." : FW Download Successful\n", "green")
        else
            WriteToLog("Device "..i.." : FW Download Failed\n", "red")
            return -1
        end
        
        
        -- RF Power Up
        if (0 == ar1.RfEnable_mult(dev_list[i])) then
            WriteToLog("Device "..i.." : RF Power Up Successful\n", "green")
        else
            WriteToLog("Device "..i.." : RF Power Up Failed\n", "red")
            return -1
        end            

        ---------------------------Static Configuration-------------------------
        
        -- Channel & ADC Configuration
        if (0 == ar1.ChanNAdcConfig_mult(dev_list[i], 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, cascade_mode_list[i])) then
            WriteToLog("Device "..i.." : Channel & ADC Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : Channel & ADC Configuration Failed\n", "red")
            return -2
        end
        
        -- Including this depends on the type of board being used.
        -- LDO configuration
        if (0 == ar1.RfLdoBypassConfig_mult(dev_list[i], 3)) then
            WriteToLog("Device "..i.." : LDO Bypass Successful\n", "green")
        else
            WriteToLog("Device "..i.." : LDO Bypass failed\n", "red")
            return -2
        end
        
        -- Low Power Mode Configuration
        if (0 == ar1.LPModConfig_mult(dev_list[i],0, 0)) then
            WriteToLog("Device "..i.." : Low Power Mode Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : Low Power Mode Configuration failed\n", "red")
            return -2
        end
        
        -- Miscellaneous Control Configuration
        if (0 == ar1.SetMiscConfig_mult(dev_list[i], 1, 0, 0, 0)) then
            WriteToLog("Device "..i.." : Misc Control Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : Misc Control Configuration failed\n", "red")
            return -2
        end        

        -- Edit this API to enable/disable the boot time calibration. Enabled by default.
        -- RF Init Calibration Configuration
        if (0 == ar1.RfInitCalibConfig_mult(dev_list[i], 1, 1, 1, 1, 1, 1, 1, 65537)) then
            WriteToLog("Device "..i.." : RF Init Calibration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : RF Init Calibration failed\n", "red")
            return -2
        end
        
        -- RF Init
        if (0 == ar1.RfInit_mult(dev_list[i])) then
            WriteToLog("Device "..i.." : RF Init Successful\n", "green")
        else
            WriteToLog("Device "..i.." : RF Init failed\n", "red")
            return -2
        end
                   
        RSTD.Sleep(500) -- A small buffer time between the devices.

        WriteToLog("AWR Device"..i.." Configuration Successful! \n", "blue")
    else
        WriteToLog("AWR Device"..i.." not enabled. Skipping! \n", "red")            
    end
    ---------------------------------------------------------------------------------
end

------------------------------------------------------------------------------
for i=1,table.getn(RadarDevice) do 

    local status    =    0        
    ----------------------Configuration of AWR Devices---------------------------
    
    if ((RadarDevice[1]==1) and (RadarDevice[i]==1)) then
    
        WriteToLog("Device"..i.." Configuration Started... \n", "blue")
        
        ---------------------------Data Configuration----------------------------------
        
        -- Data path Configuration
        if (0 == ar1.DataPathConfig_mult(dev_list[i], 0, 1, 0)) then
            WriteToLog("Device "..i.." : Data Path Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : Data Path Configuration failed\n", "red")
            return -3
        end
        
        -- Clock Configuration
        if (0 == ar1.LvdsClkConfig_mult(dev_list[i], 1, 1)) then
            WriteToLog("Device "..i.." : Clock Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : Clock Configuration failed\n", "red")
            return -3
        end
        
        -- CSI2 Configuration
        if (0 == ar1.CSI2LaneConfig_mult(dev_list[i], 1, 0, 2, 0, 4, 0, 5, 0, 3, 0, 0)) then
            WriteToLog("Device "..i.." : CSI2 Configuration Successful\n", "green")
        else
            WriteToLog("Device "..i.." : CSI2 Configuration failed\n", "red")
            return -3
        end
                   
        RSTD.Sleep(500) -- A small buffer time between the devices.

        WriteToLog("AWR Device"..i.." Configuration Successful! \n", "blue")
    else
        WriteToLog("AWR Device"..i.." not enabled. Skipping! \n", "red")            
    end
    --------------------------------------------------------------------------------
end