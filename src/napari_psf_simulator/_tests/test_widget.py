from .._widget import Psf_widget

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_psf_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()
    
    # create our widget, passing in the viewer
    my_widget = Psf_widget(viewer)
    
    assert type(my_widget.gen.wavelength) == float
    
    my_widget.calculate_psf()
    
    #assert my_widget.gen.PSF3D.shape[2] == my_widget.Nxy.val
    