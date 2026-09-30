module PersistenceDrives
parameterized_drive(p, t) = p.scale * t * (p.duration - t)
drive(t) = t * (1 - t)
duration_drive(t, duration) = t * (duration - t)
end
