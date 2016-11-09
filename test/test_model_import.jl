import Dyno

filename = joinpath(rootdir,"examples","rbc.mod")

model_data = Dyno.modfile_parser(filename)

# print(model_data)
