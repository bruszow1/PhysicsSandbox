/*
 This class was created as a rudimentary GUI for debugging.

 Most of it is copied verbatim from:
    https://docs.wxwidgets.org/latest/overview_helloworld.html
    https://wiki.wxwidgets.org/An_image_panel
 */


#include "display_window.h"
#include "wx/wxprec.h"
#include <wx/wx.h>
#include <wx/sizer.h>
#include <wx/rawbmp.h>
#include <thread>

#include "bitmap.h"
#include "movement.h"
#include "shape.h"

#ifndef WX_PRECOMP
#	include "wx/wx.h"
#endif

#define WinMain main
wxIMPLEMENT_APP(display_window);

wxBitmap bmp(IMAGE_WIDTH, IMAGE_HEIGHT, 24);
wxNativePixelData data(bmp);
int x_rotate = 0;
int y_rotate = 0;
int z_rotate = 0;
int x_increment = 0;
int y_increment = 0;
int z_increment = 0;
double timestep = 0.0;
double timestep_increment = 1.0;
bool running = false;

const int WAIT_MS = 32;

Shape **shape_arr;
unsigned int shape_count = 2;
Pixel *bitmap;

void timestep_update();

class wxImagePanel : public wxPanel
    {
        wxBitmap image;

    public:
        wxImagePanel(wxFrame* parent, wxString file, wxBitmapType format);

        wxImagePanel(wxFrame *parent, wxBitmap bmp);

        void rightClick(wxMouseEvent &event);

        void rotateEvent();
        void keyPressed(wxKeyEvent& event);


        void paintEvent(wxPaintEvent & evt);
        void paintNow();

        void render(wxDC& dc);

        // some useful events
        /*
         void mouseMoved(wxMouseEvent& event);
         void mouseDown(wxMouseEvent& event);
         void mouseWheelMoved(wxMouseEvent& event);
         void mouseReleased(wxMouseEvent& event);
         void rightClick(wxMouseEvent& event);
         void mouseLeftWindow(wxMouseEvent& event);
         void keyPressed(wxKeyEvent& event);
         void keyReleased(wxKeyEvent& event);
         */

        DECLARE_EVENT_TABLE()
    };


BEGIN_EVENT_TABLE(wxImagePanel, wxPanel)
// some useful events
/*
 EVT_MOTION(wxImagePanel::mouseMoved)
 EVT_LEFT_DOWN(wxImagePanel::mouseDown)
 EVT_LEFT_UP(wxImagePanel::mouseReleased)
 EVT_RIGHT_DOWN(wxImagePanel::rightClick)
 EVT_LEAVE_WINDOW(wxImagePanel::mouseLeftWindow)
 EVT_KEY_DOWN(wxImagePanel::keyPressed)
 EVT_KEY_UP(wxImagePanel::keyReleased)
 EVT_MOUSEWHEEL(wxImagePanel::mouseWheelMoved)
 */

// catch paint events
EVT_PAINT(wxImagePanel::paintEvent)
EVT_KEY_UP(wxImagePanel::keyPressed)
EVT_RIGHT_DOWN(wxImagePanel::rightClick)

END_EVENT_TABLE()


// some useful events
/*
 void wxImagePanel::mouseMoved(wxMouseEvent& event) {}
 void wxImagePanel::mouseDown(wxMouseEvent& event) {}
 void wxImagePanel::mouseWheelMoved(wxMouseEvent& event) {}
 void wxImagePanel::mouseReleased(wxMouseEvent& event) {}
 void wxImagePanel::rightClick(wxMouseEvent& event) {}
 void wxImagePanel::mouseLeftWindow(wxMouseEvent& event) {}
 void wxImagePanel::keyPressed(wxKeyEvent& event) {}
 void wxImagePanel::keyReleased(wxKeyEvent& event) {}
 */

wxImagePanel::wxImagePanel(wxFrame* parent, wxString file, wxBitmapType format) : wxPanel(parent) {
    // load the file... ideally add a check to see if loading was successful
    image.LoadFile(file, format);
}

wxImagePanel::wxImagePanel(wxFrame* parent, wxBitmap bmp) :
wxPanel(parent)
{
    // load the file... ideally add a check to see if loading was successful
    image = bmp;
}

/*
 * Called by the system of by wxWidgets when the panel needs
 * to be redrawn. You can also trigger this call by
 * calling Refresh()/Update().
 */

void wxImagePanel::rightClick(wxMouseEvent& event) {
    rotateEvent();
}

void wxImagePanel::keyPressed(wxKeyEvent& event) {
    if (event.GetKeyCode() == 13) {
        if (running) {
            running = false;
        } else {
            running = true;
            std::thread time_thread(timestep_update);
            time_thread.detach();
        }
    }
}

void wxImagePanel::rotateEvent()
{
    model_movement(shape_arr, shape_count, x_rotate, y_rotate, z_rotate, 0, timestep_increment);
    draw_shapes(bitmap, shape_arr, shape_count);

    wxNativePixelData::Iterator pixel(data);
    for (int i = 0; i < IMAGE_HEIGHT; i++) {
        // wxPixelData starts at upper left corner; adjust to match lower left convention
        unsigned int bitmap_y = IMAGE_HEIGHT - i - 1;
        for (int j = 0; j < IMAGE_WIDTH; j++) {
            unsigned int bitmap_index = bitmap_y * IMAGE_WIDTH + j;
            pixel.Red() = bitmap[bitmap_index].red;
            pixel.Green() = bitmap[bitmap_index].green;
            pixel.Blue() = bitmap[bitmap_index].blue;
            ++pixel;
        }
        pixel.MoveTo(data, 0, i + 1);
    }

    paintNow();
}

void wxImagePanel::paintEvent(wxPaintEvent & evt)
{
    // depending on your system you may need to look at double-buffered dcs
    wxPaintDC dc(this);
    render(dc);
}

/*
 * Alternatively, you can use a clientDC to paint on the panel
 * at any time. Using this generally does not free you from
 * catching paint events, since it is possible that e.g. the window
 * manager throws away your drawing when the window comes to the
 * background, and expects you will redraw it when the window comes
 * back (by sending a paint event).
 */
void wxImagePanel::paintNow()
{
    // depending on your system you may need to look at double-buffered dcs
    wxClientDC dc(this);
    render(dc);
}

/*
 * Here we do the actual rendering. I put it in a separate
 * method so that it can work no matter what type of DC
 * (e.g. wxPaintDC or wxClientDC) is used.
 */
void wxImagePanel::render(wxDC&  dc)
{
    dc.DrawBitmap( image, 0, 0, false );
}

enum
{
    ID_Display = 1
};

wxFrame *frame;
wxImagePanel * drawPane;

void timestep_update()
{
    while (running) {
        Sleep(WAIT_MS);

        model_movement(shape_arr, shape_count, x_rotate, y_rotate, z_rotate, 0, timestep_increment);
        draw_shapes(bitmap, shape_arr, shape_count);

        wxNativePixelData::Iterator pixel(data);
        for (int i = 0; i < IMAGE_HEIGHT; i++) {
            // wxPixelData starts at upper left corner; adjust to match lower left convention
            unsigned int bitmap_y = IMAGE_HEIGHT - i - 1;
            for (int j = 0; j < IMAGE_WIDTH; j++) {
                unsigned int bitmap_index = bitmap_y * IMAGE_WIDTH + j;
                pixel.Red() = bitmap[bitmap_index].red;
                pixel.Green() = bitmap[bitmap_index].green;
                pixel.Blue() = bitmap[bitmap_index].blue;
                ++pixel;
            }
            pixel.MoveTo(data, 0, i + 1);
        }

        drawPane->paintNow();
    }
}

// This is executed upon startup, like 'main()' in non-wxWidgets programs.
bool display_window::OnInit()
{
    // make sure to call this first
    wxInitAllImageHandlers();

    wxBoxSizer* sizer = new wxBoxSizer(wxHORIZONTAL);
    frame = new wxFrame(NULL, wxID_ANY, wxT("Display Window"), wxPoint(50,50), wxSize(IMAGE_WIDTH,IMAGE_HEIGHT));

    shape_arr = new Shape*[shape_count];
    bitmap = new Pixel[IMAGE_WIDTH * IMAGE_HEIGHT];

    shape_arr[0] = new Shape(12, 12, 9);
    double origin1[] = {200, 250, 500};
    double velocity1[] = {5, 0, 0};
    double angular_velocity1[] = {0.05, -.00, 0};
    assemble_cube45(shape_arr[0], origin1, 300, velocity1, angular_velocity1);

    shape_arr[1] = new Shape(12, 12, 9);
    double origin2[] = {650, 300, 500};
    double velocity2[] = {0, 0, 0};
    double angular_velocity2[] = {0, 0, 0};
    assemble_cube(shape_arr[1], origin2, 300, velocity2, angular_velocity2);

    model_movement(shape_arr, shape_count, 0, 0, 0, 0, timestep_increment);
    draw_shapes(bitmap, shape_arr, shape_count);

    wxNativePixelData::Iterator pixel(data);
    for (int i = 0; i < IMAGE_HEIGHT; i++) {
        // wxPixelData starts at upper left corner; adjust to match lower left convention
        unsigned int bitmap_y = IMAGE_HEIGHT - i - 1;
        for (int j = 0; j < IMAGE_WIDTH; j++) {
            unsigned int bitmap_index = bitmap_y * IMAGE_WIDTH + j;
            pixel.Red() = bitmap[bitmap_index].red;
            pixel.Green() = bitmap[bitmap_index].green;
            pixel.Blue() = bitmap[bitmap_index].blue;
            ++pixel;
        }
        pixel.MoveTo(data, 0, i + 1);
    }

    drawPane = new wxImagePanel( frame, bmp);

    sizer->Add(drawPane, 1, wxEXPAND);

    frame->SetSizer(sizer);

    frame->Show();
    return true;
}

MyFrame::MyFrame()
    : wxFrame(nullptr, wxID_ANY, "Hello World")
{
    wxMenu *menuFile = new wxMenu;
    menuFile->Append(ID_Display, "&Hello...\tCtrl-H",
                     "Help string shown in status bar for this menu item");
    menuFile->AppendSeparator();
    menuFile->Append(wxID_EXIT);

    wxMenu *menuHelp = new wxMenu;
    menuHelp->Append(wxID_ABOUT);

    wxMenuBar *menuBar = new wxMenuBar;
    menuBar->Append(menuFile, "&File");
    menuBar->Append(menuHelp, "&Help");

    SetMenuBar( menuBar );

    CreateStatusBar();
    SetStatusText("Welcome to wxWidgets!");

    Bind(wxEVT_MENU, &MyFrame::OnHello, this, ID_Display);
    Bind(wxEVT_MENU, &MyFrame::OnAbout, this, wxID_ABOUT);
    Bind(wxEVT_MENU, &MyFrame::OnExit, this, wxID_EXIT);
}

void MyFrame::OnExit(wxCommandEvent& event)
{
    running = false;

    for (unsigned int i = 0; i < shape_count; i++) {
        delete shape_arr[i];
    }
    delete [] shape_arr;
    delete [] bitmap;

    Close(true);
}

void MyFrame::OnAbout(wxCommandEvent& event)
{
    wxMessageBox("This is a wxWidgets Hello World example",
                 "About Hello World", wxOK | wxICON_INFORMATION);
}

void MyFrame::OnHello(wxCommandEvent& event)
{
    wxLogMessage("Hello world from wxWidgets!");
}