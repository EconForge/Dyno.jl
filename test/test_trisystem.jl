# import Dyno
#
# sys_str = Dict{String,String}(
#     "a"=>"1.0",
#     "b"=>"a+1",
#     "c"=>"a+b",
# )
#
# # res = Dyno.string_dict_to_symbol_dict(sys_str)
# res = Dyno.triangular_system(sys_str)
# #
# # sys = Dict(
# #     :a=>:(1.0),
# #     :b=>:(a+1),
# #     :c=>:(a+b),
# # )
