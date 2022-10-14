import numpy as np
import pandas as pd
from scipy import signal



def transform_mv_to_rcs(data_mv, amp_gain):
    """Converts data in units of mV to the internal unit representation used on 
    the RC+S device.
    
    Parameters
    ----------
    data_mv : (num_samples, num_channels) array, or transpose
        Data, either time-domain or FFT amplitude, given in units of mV. The 
        result will be returned in the same shape.
    amp_gain : int, 
        Parameter indicating the channel gain represented by the
        cfg_config_data.dev.HT_sns_ampX_gain250_trim value in the
        DeviceSettings.json file, or the metaData.ampGains OpenMind output.
    
    Returns
    -------
    data_rcs : (num_samples, num_channels) array, or transpose
        Data, either time-domain or FFT amplitude, given in internal RC+S units.
        Returned in the same shape as the input data.
    """
    
    rcs_constant = 48644.8683623726 # RCS specific conversion constant
    amp_gain = 250*(amp_gain/255);  # Convert actual channel amp gain
    data_rcs = data_mv * (amp_gain*rcs_constant) / (1000*1.2)
    
    return data_rcs


def transform_rcs_to_mv(data_rcs, amp_gain):
    """Converts data from the internal unit representation used on the RC+S 
    device to units of mV.
    
    Parameters
    ----------
    data_rcs : (num_samples, num_channels) array, or transpose
        Data, either time-domain or FFT amplitude, given in internal RC+S units.
        The result will be returned in the same shape.
    amp_gain : int, 
        Parameter indicating the channel gain represented by the
        cfg_config_data.dev.HT_sns_ampX_gain250_trim value in the
        DeviceSettings.json file, or the metaData.ampGains OpenMind output.
    
    Returns
    -------
    data_mv : (num_samples, num_channels) array, or transpose
        Data, either time-domain or FFT amplitude, given in units of mV. 
        Returned in the same shape as the input data.
    """
    
    rcs_constant = 48644.8683623726 # RCS specific conversion constant
    amp_gain = 250*(amp_gain/255);  # Convert actual channel amp gain
    data_mv = data_rcs * (1000*1.2) / (amp_gain*rcs_constant)
    
    return data_mv


def transform_ramp_rate_int_to_mas(ramp_rate_int):
    """Converts ramp rate from the integer units programmed into adaptive config
    files to units of mA/sec.
    
    Parameters
    ----------
    ramp_rate_int : positive integer
        Ramp rate, given in integer units that are used in config files.
    
    Returns
    -------
    ramp_rate_mas : float
        Ramp rate, given in units of mA/sec.
    """
    
    one_mas = 65536 # int value for 1mA/sec
    ramp_rate_mas = ramp_rate_int / one_mas
    
    return ramp_rate_mas


def transform_ramp_rate_mas_to_int(ramp_rate_mas):
    """Converts ramp rate from units of mA/sec to the integer units programmed
    into adaptive config files.
    
    Parameters
    ----------
    ramp_rate_mas : float
        Ramp rate, given in units of mA/sec.
    
    Returns
    -------
    ramp_rate_int : positive integer
        Ramp rate, given in integer units that are used in config files.
    """
    
    one_mas = 65536 # int value for 1mA/sec
    ramp_rate_int = np.around(ramp_rate_mas * one_mas).astype(int)
    
    return ramp_rate_int


def create_hann_window(L, percent=100):
    """Creates an L-point symmetric Hann window with optional modification of 
    the window peakedness.
    
    Parameters
    ----------
    L : int
        Length of the window, in units of samples
    percent : float, in the interval (0,100]. default=100
        Determines the peakedness of the window coefficients. 100% is an
        unmodified Hann window - a cosine curve with a single peak at the window
        center. Values below 100% flatten the weights by making the outer edges
        steeper and assigning a constant weight of 1 across a wider range
        centered about the window center (i.e., a plateau). As `percent` 
        approaches 0, the window approaches uniform weights.
    
    Returns
    -------
    hann_win : (N,) ndarray
        Hann window
    """
    
    if percent<=0 | percent>100:
        raise ValueError('`percent` must be in the interval (0,100]')
        
    # The actual FFT uses a smaller number of true time-domain samples and
    # zero-pads the remainder
    L = int(L)
    if L == 64:
        L_non_zero = 63 # on the device this alternates between 62 and 63
    elif L == 256:
        L_non_zero = 250
    elif L == 1024:
        L_non_zero = 1000
    else:
        L_non_zero = L
    
    # calculate the cosine curve
    B = (200/percent)*np.pi
    hann_win = 0.5*(1-np.cos(B*np.arange(L_non_zero)/(L_non_zero)))
    
    # replace the inner portion of the cosine curve with flat 1's
    if percent<100:
        peak_idx = signal.find_peaks(hann_win)[0]
        hann_win[peak_idx[0]+1:peak_idx[-1]] = 1
        
    # zero-pad remaining points
    hann_win = np.concatenate([hann_win, np.zeros(L-L_non_zero)])
        
    return hann_win
    
    
def td_to_fft(data_td, time_td, fs_td, L, interval, hann_win, 
              interp_drops=True, output_in_mv=False):
    """Computes short-time FFT of Time-Domain data as it is performed onboard 
    the RC+S device. The result may optionally be converted to match the logged 
    FFT outputs in units of mV.

    Parameters
    ----------
    data_td : (num_samples,) array
        Time-Domain data, given in internal RC+S units.
    time_td : (num_samples,) array
        Unix timestamps (in whole ms format) for the corresponding Time-Domain 
        data samples.
    fs_td : int, [250, 500, 1000]
        The Time-Domain sampling rate, in Hz.
    L : int, [64, 256, 1024]
        FFT size, in number of samples. Throws a warning if given as a value 
        that is not an option for the RC+S.
    interval : int 
        The interval, in ms, that the FFT window is shifted for producing each 
        subsequent output sample.
    hann_win : (L,) array, or transpose
        The coefficients of the Hann window applied to Time-Domain data prior to 
        computing the FFT.
    interp_drops : boolean. default=True
        Boolean flag indicating whether to linearly interpolate over missing
        Time-Domain samples. This can cause unexpected power outputs,
        particularly increases in the low frequencies, but can be advantageous 
        for downstream LD operations as it will not miss and Power Band samples.
    output_in_mv : boolean. default=False
        Boolean flag indicating whether to match the FFT output units to what is 
        logged by the device (scaled mV).

    Returns
    -------
    data_fft : (num_windows, L) array
        FFT amplitude data given in internal RC+S units, or converted to match 
        the scaled mV units that the device outputs in data logs if specified by 
        the `output_in_mv` parameter.
    time_fft : (num_windows,) array
        Unix timestamps for the corresponding FFT data samples.
    """
    
    # Check for appropriate parameter values
    if fs_td not in [250, 500, 1000]:
        raise ValueError('Time-domain sampling rate `fs_td` must be' \
                         + '250, 500, or 1000')
    if L not in [64, 256, 1024]:
        raise ValueError('FFT size `L` must be 64, 256, or 1024')
    
    # Make sure all vectors are 1-D to avoid shaping problems
    data_td = data_td.flatten()
    time_td = time_td.flatten()
    hann_win = hann_win.flatten()

    # The actual FFT uses a smaller number of true time-domain samples and
    # zero-pads the remainder
    L = int(L)
    if L == 64:
        L_non_zero = 63 # on the device this alternates between 62 and 63
    elif L == 256:
        L_non_zero = 250
    elif L == 1024:
        L_non_zero = 1000
    else:
        L_non_zero = L

    # Linearly interpolate over NaN-values
    if interp_drops:
        nan_mask = np.isnan(data_td)
        idx = np.arange(len(data_td))
        data_td[nan_mask] = np.interp(idx[nan_mask], \
                                      idx[~nan_mask], data_td[~nan_mask])

    # Pre-select all FFT window edges
    mean_window_shift = interval*fs_td/1000
    window_stops = np.uint(np.arange(L_non_zero, len(data_td), 
                                     mean_window_shift))
    window_starts = window_stops - L_non_zero
    num_windows = len(window_stops)
    time_fft = time_td[window_stops]
    data_fft = np.zeros([num_windows,L])
    # Iterate over FFT windows
    for s in range(num_windows):
        # Select the time-domain window and zero-pad remaining points
        td_window = np.zeros(L)
        td_window[:L_non_zero] = \
            data_td[window_starts[s]:window_stops[s]]
        # Apply the hann window, remove any missing datapoints, and zero-pad
        td_hann = td_window*hann_win
        td_hann = td_hann[~np.isnan(td_hann)]
        td_hann = np.concatenate([td_hann, np.zeros(L-np.shape(td_hann)[0])])
        # Take the FFT and calculate complex magnitudes
        current_fft = np.fft.fft(td_hann, L)
        current_fft = np.abs(current_fft)
        data_fft[s,:] = current_fft

    # Convert the units to match the logged output, if desired
    if output_in_mv:
        data_fft = transform_mv_to_rcs(4*data_fft/L)
    
    return data_fft, time_fft


def fft_to_pb(data_fft, fs_td, L, bit_shift, band_edges_hz=[], 
              input_is_mv=False):
    """Converts short-time FFT outputs to scaled Power Band signals (or full 
    spectrogram) with the same scaling operations performed onboard the RC+S 
    device.

    Parameters
    ----------
    data_fft : (num_windows, L) array
        FFT amplitude data given in internal RC+S units. This may also be given 
        in units of mV (matching the format in FFT data logs) if specified by 
        the `input_is_mv` parameter.
    fs_td : int, [250, 500, 1000]
        The Time-Domain sampling rate, in Hz.
    L : int, {64, 256, 1024}
        FFT size, in number of samples.
    bit_shift : int, 0:7
        Parameter indicating the number of most-significant-bits to be
        discarded. This value should be input as exactly the same value
        programmed on the device.
    band_edges_hz : optional (num_bands, 2) array, default=[]
        Edges of each power band requested, in Hz. If empty, the function will 
        return the full L/2-dimensional single-sided spectrogram.
    input_is_mv : optional boolean, default=False
        Boolean flag indicating whether the FFT input was given in units of 
        scaled mV, matching the format in the raw data logs.

    Returns
    -------
    data_pb : (num_windows, num_bands) array
        Power Band data given in internal RC+S units, or the full
        L/2-dimensional spectrogram. Note that the first bin (DC) may be
        double what it should be.
    """
    
    # Check for appropriate parameter values
    if fs_td not in [250, 500, 1000]:
        raise ValueError('Time-domain sampling rate `fs_td` must be' \
                         + '250, 500, or 1000')
    if L not in [64, 256, 1024]:
        raise ValueError('FFT size `L` must be 64, 256, or 1024')
    if bit_shift not in np.arange(8):
        raise ValueError('Bit shift parameter `bit_shift` must be in 0:7')

    # If `data_fft` was given in mV, convert back to internal RCS units
    if input_is_mv:
        data_fft = transform_mv_to_rcs(L*data_fft/4)

    # Convert amplitude to single-sided power spectrum
    data_fft = data_fft**2
    data_fft = 64 * data_fft[:,:int(L/2)] / (L**2) # all scaling steps combined
    # Perform the bit-shift
    data_pb = np.floor(data_fft/(2**(8-bit_shift)))

    # Sum over bins in each power band or return the full spectrum if none given
    if np.size(band_edges_hz)>0:
        band_edges_hz = np.reshape(band_edges_hz, (-1,2)) # Enforce shape (n,2)
        # Create a vector containing the center frequencies of all FFT bins
        center_freqs = np.arange(L/2) * fs_td/L
        # For each requested band, sum over the appropriate FFT bins
        data_pb_binned = np.zeros([np.shape(data_pb)[0], 
                                   np.shape(band_edges_hz)[0]])
        for band_idx, [min_freq, max_freq] in enumerate(band_edges_hz):
            bin_mask = (center_freqs>=min_freq) & (center_freqs<=max_freq)
            data_pb_binned[:,band_idx] = np.sum(data_pb[:,bin_mask],1)
        data_pb = data_pb_binned
    
    return data_pb


def pb_to_ld(data_pb, time_pb, update_rate, weights,
             subtract_vec=[np.zeros(4), np.zeros(4)],
             multiply_vec=[np.ones(4), np.ones(4)]):
    """Computes LD outputs from PB signals and determines state transitions.
    
    A NEW VERSION OF THIS COULD INCLUDE THE FRACTIONAL FIXED POINT VALUE WHICH
    DETERMINES THE FIXED POINT PRECISION. NOT LIKELY TO BE SUPER IMPORTANT IN
    TERMS OF ALGORITHMIC PERFORMANCE, BUT DOES SCALE THE LOGGED OUTPUTS.

    Parameters
    ----------
    data_pb : (num_pb_samples, num_bands) array
        Power Band data given in internal RC+S units.
    time_pb : (num_pb_samples, ) array
        Timestamps associated with each power band sample.
    update_rate : (2,) list of positive integers
        The number of PB samples to average over before producing an LD update.
    weights : (2,) list of (num_bands,) arrays
        LD weights. For this and all other parameters, the input may be given as
        a (2,) list to indicate the use of two LD's and their unique parameters.
    subtract_vec : optional (2,) list of (num_bands,) arrays
        default=[np.zeros(4), np.zeros(4)]
        Values directly subtracted from the PB signals prior to LD calculation.
    multiply_vec : optional (2,) list of (num_bands,) arrays
        default=[np.ones(4). np.ones(4)]
        Scaling factors applied to the PB signals after subtraction but before
        to LD calculation.

    Returns
    -------
    ld_output : (2,) list of (num_ld_updates,) arrays
        Continuous-valued LD outputs.
    time_ld : (2,) list of (num_ld_updates,) arrays
        Timestamps associated with each LD update sample.
    update_tbl : (num_total_ld_updates, 3) array
        A sorted lookup table for indexing outputs from each LD in order during 
        stream processing. Updates from both LD's are represented in the same 
        array. First column is the PB sample index (shared clock for both LD's),
        second column is the LD sample index (unique to each LD), and third
        column is the LD identifier (0/1).
    """
    
    num_pb_samples, num_bands = np.shape(data_pb)
    if np.size(weights[1]) == 0:
        num_lds = 1
    else:
        num_lds = 2
    
    # Test for improper inputs and reformat inputs if necessary
    for k in range(num_lds):
        if update_rate[k]==0:
            update_rate[k]=1
    
    # Compute LD outputs and prepare for stream processing
    ld_output = [[],[]]
    time_ld = [[],[]]
    update_tbl = np.empty([0,3])
    for k in range(num_lds):
        # Normalize the PB signals
        subtract_mat = np.tile(subtract_vec[k][:num_bands], (num_pb_samples,1))
        multiply_mat = np.tile(multiply_vec[k][:num_bands], (num_pb_samples,1))
        pb_norm = np.multiply(data_pb-subtract_mat, multiply_mat)
        
        # Compute LD outputs without blanking. Blanking will be done in stream.
        clip_samples = int(num_pb_samples / update_rate[k]) * update_rate[k]
        pb_norm = pb_norm[:clip_samples, :]
        pb_reshaped = np.reshape(pb_norm, [-1, update_rate[k], num_bands])
        pb_updated = np.mean(pb_reshaped, axis=1)
        ld_output[k] = np.dot(pb_updated, np.reshape(weights[k], [-1,1]))
        
        # Assign timestamp array for the LD outputs
        pb_idx = np.arange(update_rate[k]-1, clip_samples, update_rate[k])
        time_ld[k] = time_pb[pb_idx]
        
        # Create a sorted lookup table for indexing outputs from each LD
        ld_sample_idx = np.arange(np.shape(ld_output[k])[0])
        ld_id = k*np.ones(np.shape(ld_output[k]))
        single_ld_updates = np.concatenate((pb_idx[:,np.newaxis], 
                                            ld_sample_idx[:,np.newaxis],
                                            ld_id), axis=1)
        update_tbl = np.append(update_tbl, single_ld_updates, axis=0)
    update_tbl = update_tbl[update_tbl[:,0].argsort(),:].astype(int)
   
    return ld_output, time_ld, update_tbl


def ld_to_state(ld_output, update_tbl, time_pb, update_rate, dual_threshold, 
                threshold, onset_duration, termination_duration, blank_duration, 
                blank_both=[False, False]):
    """Determines LD state transitions from the LD outputs.

    Parameters
    ----------
    ld_output : (2,) list of (num_ld_updates,) arrays
        Continuous-valued LD outputs.
    update_tbl : (num_total_ld_updates, 3) array
        A sorted lookup table for indexing outputs from each LD in order during 
        stream processing. Updates from both LD's are represented in the same 
        array. First column is the PB sample index (shared clock for both LD's),
        second column is the LD sample index (unique to each LD), and third
        column is the LD identifier (0/1).
    time_pb : (num_pb_samples, ) array
        Timestamps associated with each power band sample.
    update_rate : (2,) list of positive integers
        The number of PB samples to average over before producing an LD update.
    weights : (2,) list of (num_bands,) arrays
        LD weights. For this and all other parameters, the input may be given as
        a (2,) list to indicate the use of two LD's and their unique parameters.
    dual_threshold : (2,) list of booleans
        Indicates whether each LD is using two thresholds or one. If a single
        threshold is being used, LD states will be returned as 0 or 1.
    threshold : (2,) list of {integers or (2,) arrays}
        Lower and upper thresholds for each LD. If dual_threshold is True, two
        thresholds must be given.
    onset_duration : (2,) list of positive integers
        The number of LD updates (or outputs) that must be above the threshold
        in order to trigger a state change to the super-threshold state.
    termination_duration : (2,) list of positive integers
        The number of LD updates (or outputs) that must be below the threshold
        in order to trigger a state change to the sub-threshold state.
    blank_duration : (2,) list of positive integers
        The number of PB samples following a state change that each LD will hold
        off changing state again.
    blank_both : optional (2,) boolean array; default=[False, False]
        Indicates whether both LD's will be blanked when a single LD triggers
        state change blanking. The duration of the blanking will be determined
        by the state_blank corresponding to the LD that triggered blanking.

    Returns
    -------
    state : (num_pb_samples/update_rate,) array
        Discrete LD state at each timepoint. States can be 0:8.
    time_state : (num_ld_updates,) array
        Timestamps associated with each LD update sample.
    ld_output : (2,) list of (num_ld_updates,) arrays
        Continuous-valued LD outputs, updated to be "frozen" during blanking.
    """
    
    if np.size(ld_output[1]) == 0:
        num_lds = 1
    else:
        num_lds = 2
    
    # Test for improper inputs and reformat inputs if necessary
    for k in range(num_lds):
        if onset_duration[k]==0:
            onset_duration[k]=1
        if termination_duration[k]==0:
            termination_duration[k]=1
        if blank_duration[k]==0:
            blank_duration[k]=1
        
    # Update LD states using stream processing
    state_history = [np.zeros(np.max([onset_duration[k],
                                      termination_duration[k]]))
                     for k in range(num_lds)]
    blank_counter = [0,0]
    state = np.empty([0]).astype(int)
    current_state = [0,0]
    for idx, update in enumerate(update_tbl):
        # Log info about current update
        ld_sample_idx = update[1]
        k = update[2]
        cur_ld_output = ld_output[k][ld_sample_idx]
        
        # Repeat LD output and state if currently blanked
        if blank_counter[k]>0: 
            ld_output[k][ld_sample_idx] = ld_output[k][ld_sample_idx-1]
            state = np.append(state, state[-1]) 
            fresh_blank = False
        
        # Compute new state if not blanked
        else:       
            # Update the recent history of "immediate" ld states
            state_history[k] = np.roll(state_history[k], 1)
            if dual_threshold[k]:
                state_history[k][0] = int(cur_ld_output>threshold[k][0]) \
                                      + int(cur_ld_output>threshold[k][1])
            else:
                state_history[k][0] = int(cur_ld_output>threshold[k])

            # Update the single LD state
            current_state[k], blank_counter[k], fresh_blank = \
                                determine_single_ld_current_state(
                                     state_history[k], current_state[k], 
                                     onset_duration[k], termination_duration[k], 
                                     blank_duration[k], blank_counter[k]) 
            if fresh_blank and blank_both[k] and num_lds==2:
                blank_counter[np.abs(1-k)] = blank_duration[np.abs(1-k)]

            # Update the combined LD state
            state = np.append(state, current_state[0]+3*current_state[1])
        
        # Update the blanking counter
        if idx+1 >= np.shape(update_tbl)[0]:
            d_interval = 0
        elif fresh_blank:
            d_interval = 1
        else:
            d_interval = update_tbl[idx+1,0] - update_tbl[idx,0]
        blank_counter[0] -= d_interval
        if num_lds==2:
            blank_counter[1] -= d_interval
        
    # Create the timestamp vector
    time_state = time_pb[update_tbl[:,0]]
    
    # Correct for situation where both LD's updated on the same sample,
    # resulting in two state changes at the same moment in time. Return only
    # the combined LD state AFTER both of the LD's had updated.
    time_state, idx = np.unique(np.flip(time_state), return_index=True)
    state = np.flip(state)[idx]
   
    return state, time_state, ld_output


def determine_single_ld_current_state(state_history, prev_state, 
                                      onset_duration, termination_duration, 
                                      blank_duration, blank_counter):
    """Applies the onset/termination duration to determine the current LD state.

    Parameters
    ----------
    state_history : (num_ld_samples,) array
        The recent history of "immediate" LD states, indicating which side of
        the thresholds the LD output was on at each instant.
    prev_state : positive integer, 0:8
        The previous confirmed (as opposed to immediate) state of the LD.
    onset_duration : positive integer
        The number of LD updates (or outputs) that must be above the threshold
        in order to trigger a state change to the super-threshold state.
    termination_duration : positive integer
        The number of LD updates (or outputs) that must be below the threshold
        in order to trigger a state change to the sub-threshold state.
    blank_duration : positive integer
        The number of PB samples following a state change that each LD will hold
        off changing state again.
    blank_counter : integer
        The number of remaining PB samples that the LD state change will be 
        blanked.

    Returns
    -------
    prev_state : positive integer, 0:8
        The current confirmed (as opposed to immediate) state of the LD.
    blank_counter : integer
        The number of remaining PB samples that the LD state change will be 
        blanked.
    fresh_blank : boolean
        Shows whether the LD has just changed states and initiated blanking
    """
    
    fresh_blank = False
    if np.all(state_history[:onset_duration]==2): # {0,1} -> 2
        current_state = 2
    elif np.all(state_history[:termination_duration]==0): # {1,2} -> 0
        current_state = 0
    elif (prev_state==2) \
         & np.all(state_history[:termination_duration]<2): # 2 -> 1
        current_state = 1
    elif (prev_state==0) \
         & np.all(state_history[:onset_duration]>0): # 0 -> 1
        current_state = 1
    else:
        current_state = prev_state
        
    if current_state != prev_state:
        blank_counter = blank_duration
        fresh_blank = True
        
    return current_state, blank_counter, fresh_blank


def state_to_stim(state, time_state, target_amp, rise_time, fall_time):
    """Predicts stimulation amplitude time series from LD states, target
    amplitudes, and rise/fall times.

    Parameters
    ----------
    state : (num_ld_updates,) array
        Discrete LD state at each timepoint. States can be 0:8.
    time_state : (num_ld_updates,) array
        Timestamps associated with each LD update sample.
    target_amp : (8,) array
        The target stimulation amplitude for each state. A value of 25.5 will
        cause the stimulation amplitude to remain at the same value it was upon
        state change.
    rise_time : positive integer
        Increasing stim ramp rate, given in units of mA/sec.
    fall_time : positive integer
        Decreasing stim ramp rate, given in units of mA/sec.

    Returns
    -------
    stim : (num_stim_updates,) array
        The stimulation amplitude, in mA. Only logged at change points 
        (piecewise linear function)
    time_stim : (num_stim_updates,) array
        The timestamp for each stim sample.
    """
    
    state_change_idx = np.squeeze(np.argwhere(np.diff(state)!=0) + 1)
    stim = np.array([0, target_amp[state[0]]])
    time_stim = np.array([time_state[0], time_state[0] + stim[-1] / rise_time])
    for idx in state_change_idx:
        if time_stim[-1] > time_state[idx]: #if state changes during a stim ramp
            # replace the last forecasted sample with an interpolated one
            time_stim[-1] = time_state[idx]
            if stim[-1] > stim[-2]: # was ramping up
                stim[-1] = stim[-2] \
                           + (time_stim[-1] - time_stim[-2]) * rise_time
            else: # was ramping down
                stim[-1] = stim[-2] \
                           - (time_stim[-1] - time_stim[-2]) * fall_time
            # forecast the end of the ramp
            if target_amp[state[idx]] < 25.5: # if not holding, then begin ramp
                stim = np.append(stim, target_amp[state[idx]])
                if stim[-1] > stim[-2]:
                    ramp_duration = (stim[-1] - stim[-2]) / rise_time
                else:
                    ramp_duration = (stim[-2] - stim[-1]) / fall_time
                time_stim = np.append(time_stim, time_stim[-1] + ramp_duration)
        else: # if state changes during steady-state stim
            # report the preceding stim amp
            time_stim = np.append(time_stim, time_state[idx])
            stim = np.append(stim, stim[-1])
            # forecast the end of the ramp
            if target_amp[state[idx]] < 25.5: # if not holding, then begin ramp
                stim = np.append(stim, target_amp[state[idx]])
                if stim[-1] > stim[-2]:
                    ramp_duration = (stim[-1] - stim[-2]) / rise_time
                else:
                    ramp_duration = (stim[-2] - stim[-1]) / fall_time
                time_stim = np.append(time_stim, time_stim[-1] + ramp_duration)
    if time_stim[-1] < time_state[-1]: # add a final endpoint
        time_stim = np.append(time_stim, time_state[-1])
        stim = np.append(stim, stim[-1])
    
    return stim, time_stim
