from napari_psf_simulator import Psf_widget

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_PSF_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()
    
    # create our widget, passing in the viewer
    my_widget = Psf_widget(viewer)

    # call our widget method
    #my_widget.calculate_psf()

    # read captured output and check that it's as we expected
    #captured = capsys.readouterr()
    #assert captured.out == "napari has 1 layers\n"
    