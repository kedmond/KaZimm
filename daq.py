# This is a near-verbatim translation of the example program found at 
# http://www.scipy.org/Cookbook/Data_Acquisition_with_NIDAQmx
import ctypes
import numpy
def getData(sampRate, duration):
    nidaq = ctypes.windll.nicaiu # load the DLL
##############################
    # Setup some typedefs and constants
    # to correspond with values in
    # C:\Program Files\National Instruments\NI-DAQ\DAQmx ANSI C Dev\include\NIDAQmx.h
    # the typedefs
    int32 = ctypes.c_long
    uInt32 = ctypes.c_ulong
    uInt64 = ctypes.c_ulonglong
    float64 = ctypes.c_double
    TaskHandle = uInt32
    # the constants
    DAQmx_Val_Cfg_Default = int32(-1)
    DAQmx_Val_Volts = 10348
    DAQmx_Val_Rising = 10280
    DAQmx_Val_FiniteSamps = 10178
    DAQmx_Val_GroupByChannel = 0
##############################
    def CHK(err):
        """a simple error checking routine"""
        if err < 0:
            buf_size = 10000
            buf = ctypes.create_string_buffer('\000' * buf_size)
            nidaq.DAQmxGetErrorString(err,ctypes.byref(buf),buf_size)
            raise RuntimeError('nidaq call failed with error %d: %s'%(err,repr(buf.value)))
    # initialize variables
    taskHandle = TaskHandle(0)
    max_num_samples = int(duration*sampRate)
    data = numpy.zeros((max_num_samples,),dtype=numpy.float64)
    # now, on with the program
    CHK(nidaq.DAQmxCreateTask("",ctypes.byref(taskHandle)))
    CHK(nidaq.DAQmxCreateAIVoltageChan(taskHandle,"Dev1/ai0","",DAQmx_Val_Cfg_Default,float64(-10.0),float64(10.0),DAQmx_Val_Volts,None))
    CHK(nidaq.DAQmxCfgSampClkTiming(taskHandle,"",float64(sampRate),DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,uInt64(max_num_samples)));
    CHK(nidaq.DAQmxStartTask(taskHandle))
    read = int32()
    CHK(nidaq.DAQmxReadAnalogF64(taskHandle,max_num_samples,float64(duration),DAQmx_Val_GroupByChannel,data.ctypes.data,max_num_samples,ctypes.byref(read),None))
# Kazem - 6/Nov/2013 - I changed float64(10.0) to float64(duration).
# 10.0 limited the acquisitoin period to 10 seconds
#  I do not know the rationale for this, but it works really well now.
#    print "Acquired %d points"%(read.value)
    nidaq.DAQmxStopTask(taskHandle)
    nidaq.DAQmxClearTask(taskHandle)
    return data