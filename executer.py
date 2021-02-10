import fit_real_refs_2layers
import winsound
frequency = 1000  # Set Frequency To 2500 Hertz
duration = 400  # Set Duration To 1000 ms == 1 second
for i in range(10):
    winsound.Beep(frequency, duration)