struct QPSolveFailure <: Exception
    terminationstatus
    primalstatus
    dualstatus
end

Base.showerror(io::IO, err::QPSolveFailure) = print(io, """
    QP solve unsuccessful.
    Termination status: $(err.terminationstatus)
    Primal status: $(err.primalstatus)
    Dual status: $(err.dualstatus)""")
