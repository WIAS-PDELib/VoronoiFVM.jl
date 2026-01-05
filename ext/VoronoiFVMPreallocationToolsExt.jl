module VoronoiFVMPreallocationToolsExt

using PreallocationTools: DiffCache
using VoronoiFVM: VoronoiFVM

VoronoiFVM.myround(c::DiffCache; kwargs...) = string(nameof(typeof(c)))

end
