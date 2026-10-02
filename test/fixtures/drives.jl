module PersistenceDrives
parameterized_drive(p, t, duration) = p.scale * t * (duration - t)
drive(t) = t * (1 - t)
duration_drive(t, duration) = t * (duration - t)
carrier(p, t) = sinpi(2p.frequency * t + p.phase / pi)
end
