# Images

All reconstructed data is stored as an [AxisArray](https://github.com/JuliaArrays/AxisArrays.jl).
The `AxisArrays` package is part of the [Images](https://juliaimages.github.io/latest/) package
family, which groups all image processing related functionality together. We note that the
term `Image` does not restrict the dimensionality of the data types to 2D but in fact images
can be of arbitrary dimensionality.

The reconstructed MRI image `I` is an `AxisArray` and has five dimensions. The first three
are the spatial dimension `x`, `y`, and `z`, whereas dimension four encodes the number
of echos that have been reconstructed, while dimension five encodes individual coils that
may have been reconstructed independently. By using an `AxisArray` the object does not
only consist of the data but it additionally encodes the physical size of the image
as well as the echo times. To extract the ordinary Julia array one can simply use `Ireco.data`.

The advantage of encoding the physical dimensions is the image data can be stored
without loosing the dimensions of the data. For instance one can call
```julia
saveImage(filename, I)
```
to store the image and
```julia
I = loadImage(filename)
```
to load the image. Currently, MRIReco does support the NIfTI file format. By default,
`saveImage` stores the data complex valued if the image `I` is complex valued.
To store the magnitude image one can call
```julia
saveImage(filename, I, true)
```
