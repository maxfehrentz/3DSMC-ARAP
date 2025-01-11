#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkCommand.h>
#include <vtkGlyph3DMapper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkOBJReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static MouseInteractorStyle* New();
    vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    MouseInteractorStyle()
    {
        SelectedVertexId = -1;
        IsDragging = false;
    }

    virtual void OnKeyPress() override
    {
        vtkRenderWindowInteractor* rwi = this->Interactor;
        std::string key = rwi->GetKeySym();
        std::string keyState = rwi->GetKeySym();  // Get the key state

        if (key == "x" || key == "X")
        {
            std::cout << "OnKeyPress - Key: " << key 
                      << " RepeatCount: " << rwi->GetRepeatCount()
                      << " KeyCode: " << rwi->GetKeyCode()
                      << std::endl;
            // Get current mouse position
            int* clickPos = this->GetInteractor()->GetEventPosition();

            vtkSmartPointer<vtkCellPicker> picker =
                vtkSmartPointer<vtkCellPicker>::New();
            picker->SetTolerance(0.0005);

            picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
            vtkActor* actor = picker->GetActor();

            if (actor == this->MeshActor)
            {
                vtkIdType id = picker->GetPointId();
                if (id != -1)
                {
                    SelectedVertexId = id;
                    double* pos = this->MeshPolyData->GetPoint(id);
                    originalPosition[0] = pos[0];
                    originalPosition[1] = pos[1];
                    originalPosition[2] = pos[2];
                    IsDragging = true;
                    this->GetInteractor()->GetRenderWindow()->SetDesiredUpdateRate(10.0);
                    return;
                }
            }
        }
        else if (key == "Escape" && IsDragging)
        {
            std::cout << "Was dragging and Escape key pressed" << std::endl;
            // Cancel vertex manipulation
            IsDragging = false;
            SelectedVertexId = -1;
        }
        
        // Don't forget to forward other keys to the parent class
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

    virtual void OnMouseMove() override
    {
        if (IsDragging && SelectedVertexId != -1)
        {
            int* mousePos = this->GetInteractor()->GetEventPosition();
            double worldPos[4];
            this->GetDefaultRenderer()->SetDisplayPoint(mousePos[0], mousePos[1], 0);
            this->GetDefaultRenderer()->DisplayToWorld();
            this->GetDefaultRenderer()->GetWorldPoint(worldPos);
            if (worldPos[3] != 0.0)
            {
                for (int i = 0; i < 3; i++)
                {
                    worldPos[i] /= worldPos[3];
                }
                double newPos[3] = { worldPos[0], worldPos[1], originalPosition[2] };
                this->MeshPolyData->GetPoints()->SetPoint(SelectedVertexId, newPos);
                this->MeshPolyData->GetPoints()->Modified();
                this->MeshMapper->Update();
                this->GetInteractor()->GetRenderWindow()->Render();
            }
            return;
        }

        // Always call parent's OnMouseMove to maintain camera navigation
        vtkInteractorStyleTrackballCamera::OnMouseMove();
    }

    void SetMeshActor(vtkActor* actor)
    {
        this->MeshActor = actor;
    }

    void SetMeshPolyData(vtkPolyData* polyData)
    {
        this->MeshPolyData = polyData;
    }

    void SetMeshMapper(vtkPolyDataMapper* mapper)
    {
        this->MeshMapper = mapper;
    }

private:
    vtkActor* MeshActor;
    vtkPolyData* MeshPolyData;
    vtkPolyDataMapper* MeshMapper;
    vtkIdType SelectedVertexId;
    double originalPosition[3];
    bool IsDragging;
};

vtkStandardNewMacro(MouseInteractorStyle);

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " Filename.obj" << std::endl;
        return EXIT_FAILURE;
    }

    // Read the OBJ file
    vtkSmartPointer<vtkOBJReader> reader =
        vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    // Create mapper and actor for the mesh
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> meshActor =
        vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(mapper);
    meshActor->GetProperty()->SetColor(0.8, 0.8, 0.8);

    // Create sphere glyphs for vertices
    vtkSmartPointer<vtkSphereSource> sphereSource =
        vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.001);
    sphereSource->Update();

    vtkSmartPointer<vtkGlyph3DMapper> glyphMapper =
        vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyphMapper->SetInputData(reader->GetOutput());
    glyphMapper->SetSourceConnection(sphereSource->GetOutputPort());
    glyphMapper->Update();

    vtkSmartPointer<vtkActor> glyphActor =
        vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(1.0, 0.0, 0.0);

    // Renderer setup
    vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(meshActor);
    renderer->AddActor(glyphActor);
    renderer->SetBackground(0.1, 0.2, 0.4);

    // Render window
    vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("Interactive OBJ Viewer");

    // Interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Set custom interactor style
    vtkSmartPointer<MouseInteractorStyle> style =
        vtkSmartPointer<MouseInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->SetMeshActor(meshActor);
    style->SetMeshPolyData(reader->GetOutput());
    style->SetMeshMapper(mapper);
    renderWindowInteractor->SetInteractorStyle(style);

    // Start interaction
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}