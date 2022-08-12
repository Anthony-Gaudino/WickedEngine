#include "stdafx.h"
#include "AnimationWindow.h"
#include "Editor.h"

using namespace wi::ecs;
using namespace wi::scene;

void AnimationWindow::Create(EditorComponent* _editor)
{
	editor = _editor;
	wi::gui::Window::Create(ICON_ANIMATION " Animation", wi::gui::Window::WindowControls::COLLAPSE);
	SetSize(XMFLOAT2(520, 140));

	float x = 80;
	float y = 0;
	float hei = 18;
	float wid = 200;
	float step = hei + 2;

	loopedCheckBox.Create("Looped: ");
	loopedCheckBox.SetTooltip("Toggle animation looping behaviour.");
	loopedCheckBox.SetSize(XMFLOAT2(hei, hei));
	loopedCheckBox.SetPos(XMFLOAT2(x, y += step));
	loopedCheckBox.OnClick([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			animation->SetLooped(args.bValue);
		}
	});
	AddWidget(&loopedCheckBox);

	playButton.Create("Play");
	playButton.SetTooltip("Play/Pause animation.");
	playButton.SetSize(XMFLOAT2(100, hei));
	playButton.SetPos(XMFLOAT2(loopedCheckBox.GetPos().x + loopedCheckBox.GetSize().x + 5, y));
	playButton.OnClick([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			if (animation->IsPlaying())
			{
				animation->Pause();
			}
			else
			{
				animation->Play();
			}
		}
	});
	AddWidget(&playButton);

	stopButton.Create("Stop");
	stopButton.SetTooltip("Stop animation.");
	stopButton.SetSize(XMFLOAT2(100, hei));
	stopButton.SetPos(XMFLOAT2(playButton.GetPos().x + playButton.GetSize().x + 5, y));
	stopButton.OnClick([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			animation->Stop();
		}
	});
	AddWidget(&stopButton);

	timerSlider.Create(0, 1, 0, 100000, "Timer: ");
	timerSlider.SetSize(XMFLOAT2(wid, hei));
	timerSlider.SetPos(XMFLOAT2(x, y += step));
	timerSlider.OnSlide([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			animation->timer = args.fValue;
		}
	});
	timerSlider.SetEnabled(false);
	timerSlider.SetTooltip("Set the animation timer by hand.");
	AddWidget(&timerSlider);

	amountSlider.Create(0, 1, 1, 100000, "Amount: ");
	amountSlider.SetSize(XMFLOAT2(wid, hei));
	amountSlider.SetPos(XMFLOAT2(x, y += step));
	amountSlider.OnSlide([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			animation->amount = args.fValue;
		}
	});
	amountSlider.SetEnabled(false);
	amountSlider.SetTooltip("Set the animation blending amount by hand.");
	AddWidget(&amountSlider);

	speedSlider.Create(0, 4, 1, 100000, "Speed: ");
	speedSlider.SetSize(XMFLOAT2(wid, hei));
	speedSlider.SetPos(XMFLOAT2(x, y += step));
	speedSlider.OnSlide([&](wi::gui::EventArgs args) {
		AnimationComponent* animation = editor->GetCurrentScene().animations.GetComponent(entity);
		if (animation != nullptr)
		{
			animation->speed = args.fValue;
		}
	});
	speedSlider.SetEnabled(false);
	speedSlider.SetTooltip("Set the animation speed.");
	AddWidget(&speedSlider);



	SetMinimized(true);
	SetVisible(false);

}

void AnimationWindow::SetEntity(Entity entity)
{
	this->entity = entity;
}

void AnimationWindow::Update()
{
	Scene& scene = editor->GetCurrentScene();

	if (!scene.animations.Contains(entity))
	{
		SetEnabled(false);
		return;
	}
	else
	{
		SetEnabled(true);
	}

	AnimationComponent& animation = *scene.animations.GetComponent(entity);

	if (animation.IsPlaying())
	{
		playButton.SetText("Pause");
	}
	else
	{
		playButton.SetText("Play");
	}

	loopedCheckBox.SetCheck(animation.IsLooped());

	timerSlider.SetRange(0, animation.GetLength());
	timerSlider.SetValue(animation.timer);
	amountSlider.SetValue(animation.amount);
	speedSlider.SetValue(animation.speed);
}
