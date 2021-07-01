
class SetVisibilityCallback:
    """Helper callback to keep a reference to the actor being modified."""

    def __init__(self, actors):
        if type(actors) is list:
            self.actors = actors
        else:
            self.actors = [actors]

    def __call__(self, state):
        for _actor in self.actors:
            if type(_actor) is vtkmodules.vtkInteractionWidgets.vtkSliderWidget:
                pass
            else:
                _actor.SetVisibility(state)