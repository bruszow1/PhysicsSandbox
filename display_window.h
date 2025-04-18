#ifndef DISPLAY_WINDOW_H
#define DISPLAY_WINDOW_H
#include <wx/wx.h>

class display_window : public wxApp
{
    public:
        virtual bool OnInit() override;

};

class MyFrame : public wxFrame
{
public:
    MyFrame();

private:
    void OnHello(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
};

#endif //DISPLAY_WINDOW_H
