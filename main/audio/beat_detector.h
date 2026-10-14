#ifndef BEAT_DETECTOR_H
#define BEAT_DETECTOR_H

#include <vector>
#include <cstdint>

/**
 * 音乐节拍检测器（律动模式用）。
 *
 * 从 finch-chan 的 AudioVisualizer 移植的简化版：
 *   - 512 点 FFT + Hann 窗，取幅度谱
 *   - 只盯低频鼓点段（约 40~160Hz 的 12 个频段），低频能量突增 = 一个鼓点
 *   - 最小节拍间隔 240ms（约 250 BPM 上限），避免一个鼓点重复触发
 *
 * 输入：16kHz 单声道 int16 PCM（每次喂一段，内部累积到 512 点）
 * 输出：检测到节拍返回 true（一次只报一个节拍）
 */
class BeatDetector {
public:
    static constexpr uint16_t kFftSize = 512;
    static constexpr uint16_t kSampleRate = 16000;

    BeatDetector();

    /**
     * 喂一段 PCM 样本（16kHz int16 单声道），内部累积到 512 点后算一次节拍。
     * 返回 true 表示这一帧里检测到一个节拍（鼓点）。
     */
    bool Process(const std::vector<int16_t>& pcm);

    /** 清空内部缓冲（切歌/停止律动时调用） */
    void Reset();

private:
    // 512 点 FFT 缓冲（累积 PCM）
    std::vector<float> buffer_;
    size_t buffer_pos_ = 0;

    // Hann 窗
    std::vector<float> window_;

    // 上一个节拍的时间（毫秒），用于最小间隔限制
    int64_t last_beat_ms_ = 0;

    // 低频能量的滑动均值（用于"突增"判定：当前低频能量远超历史均值 = 一个鼓点）
    float avg_low_energy_ = 0.0f;

    bool ComputeBeat(const std::vector<float>& samples, int64_t now_ms);
};

#endif // BEAT_DETECTOR_H
