# Same-filesystem rename commit; a failed replacement retains/restores the old directory.
function _commit_bundle!(stage,destination;overwrite=false,rename=mv)
    ispath(destination) && !overwrite && _ioerror(destination,"destination exists; use overwrite=true")
    islink(destination) && _ioerror(destination,"destination symlink is forbidden")
    ispath(destination) && !isdir(destination) && _ioerror(destination,"destination is not a directory")
    backup=nothing
    if isdir(destination)
        backup=tempname(dirname(destination))
        rename(destination,backup)
    end
    try
        rename(stage,destination)
    catch original
        if !isnothing(backup)
            try
                rename(backup,destination)
            catch rollback
                throw(CompositeException(original,ErrorException("rollback failed; original directory remains at $backup: $rollback")))
            end
        end
        rethrow()
    end
    !isnothing(backup) && rm(backup;recursive=true)
    return destination
end
